"""A small reader and writer for OpenFOAM dictionary files.

The reader covers the subset RDEsim configurations use: nested dictionaries,
lists, strings, numbers, words, ``//`` and ``/* */`` comments, and
``#include``. It is not a full Foam dictionary implementation.
"""

from __future__ import annotations

import re
from pathlib import Path


class FoamParseError(Exception):
    pass


class FoamDict(dict):
    """Ordered dictionary produced by the reader."""


class FoamList(list):
    """Parenthesised list produced by the reader."""


class FoamTokens(list):
    """Several tokens before a semicolon, written back without parentheses."""


def parse_file(path: str | Path) -> FoamDict:
    path = Path(path)
    text = path.read_text(encoding="utf-8")
    return parse(text, source=str(path), base_dir=path.parent)


def parse(text: str, source: str = "<string>", base_dir: Path | None = None) -> FoamDict:
    parser = _Parser(text, source, base_dir or Path("."))
    data = FoamDict()
    while parser.peek() is not None:
        parser.parse_entry(data)
    return data


def write_file(path: str | Path, class_name: str, object_name: str, body: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = foam_header(class_name, object_name, path.parent.name)
    path.write_text(header + body.rstrip() + "\n\n// ************************************************************************* //\n", encoding="utf-8")


def foam_header(class_name: str, object_name: str, location: str) -> str:
    return f"""/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  14
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {class_name};
    location    "{location}";
    object      {object_name};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def format_number(value: float | int) -> str:
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, int) and not isinstance(value, bool):
        return str(value)
    text = format(float(value), ".12g")
    if text in ("-0", "-0.0"):
        return "0"
    return text


def format_data(value) -> str:
    if isinstance(value, FoamDict) or isinstance(value, dict):
        return format_dict(value, indent=0)
    if isinstance(value, FoamTokens):
        return " ".join(format_data(item) for item in value)
    if isinstance(value, (FoamList, list, tuple)):
        return format_list(value)
    if isinstance(value, bool):
        return "yes" if value else "no"
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return format_number(value)
    if isinstance(value, str):
        if re.fullmatch(r"[A-Za-z_][\w:.+-]*", value):
            return value
        return '"' + value.replace('"', '\\"') + '"'
    raise TypeError(f"Cannot write value of type {type(value).__name__}")


def render_entries(data, indent: int = 0) -> str:
    pad = "    " * indent
    lines = []
    for key, value in data.items():
        if isinstance(value, (FoamDict, dict)):
            lines.append(f"{pad}{key}")
            lines.append(format_dict(value, indent))
        elif isinstance(value, (FoamList, list, tuple)) and value and isinstance(value[0], (dict, FoamDict)):
            lines.append(f"{pad}{key}")
            lines.append(format_dict_list(value, indent) + ";")
        else:
            lines.append(f"{pad}{key}    {format_data(value)};")
    return "\n".join(lines)


def format_dict(data, indent: int = 0) -> str:
    pad = "    " * indent
    body = render_entries(data, indent + 1)
    return f"{pad}{{\n{body}\n{pad}}}"


def format_dict_list(items, indent: int) -> str:
    pad = "    " * indent
    lines = [f"{pad}("]
    for item in items:
        rendered = format_dict(item, indent + 1)
        lines.extend(pad + "    " + line for line in rendered.splitlines())
    lines.append(f"{pad})")
    return "\n".join(lines)


def format_list(value) -> str:
    parts = [format_data(item) if not isinstance(item, (dict, FoamDict)) else format_dict(item) for item in value]
    if any(isinstance(item, (dict, FoamDict)) or (isinstance(item, str) and "\n" in item) for item in value):
        return "(\n" + "\n".join("    " + part for part in parts) + "\n)"
    return "(" + " ".join(parts) + ")"


class _Parser:
    def __init__(self, text: str, source: str, base_dir: Path):
        self.text = text
        self.source = source
        self.base_dir = base_dir
        self.n = len(text)
        self.i = 0

    def peek(self):
        self._skip()
        if self.i >= self.n:
            return None
        return self.text[self.i]

    def parse_entry(self, dest: FoamDict) -> None:
        self._skip()
        if self.i >= self.n:
            return
        if self.text.startswith("#include", self.i):
            self._include(dest)
            return
        key = self._word()
        self._skip()
        if self.i >= self.n:
            raise FoamParseError(f"{self.source}: key '{key}' has no value")
        value = self._value_sequence()
        dest[key] = value

    def _include(self, dest: FoamDict) -> None:
        self.i += len("#include")
        self._skip()
        if self.peek() != '"':
            raise FoamParseError(f"{self.source}: #include expects a quoted path")
        raw = self._string()
        path = Path(raw)
        if not path.is_absolute():
            path = self.base_dir / path
        if not path.is_file():
            raise FoamParseError(f"{self.source}: #include file not found: {path}")
        included = parse_file(path)
        for key, value in included.items():
            if key == "FoamFile":
                continue
            dest[key] = value
        self._optional_semicolon()

    def _value_sequence(self):
        self._skip()
        if self.peek() in "{(":
            value = self._compound()
            self._optional_semicolon()
            return value
        first = self._token()
        pieces = [first]
        while True:
            self._skip()
            if self.i >= self.n or self.text[self.i] in ";})":
                break
            if self.text[self.i] in "{(":
                pieces.append(self._compound())
            else:
                pieces.append(self._token())
        self._expect_semicolon()
        if len(pieces) == 1:
            return pieces[0]
        return FoamTokens(pieces)

    def _compound(self):
        self._skip()
        char = self.peek()
        if char == "{":
            return self._dict()
        if char == "(":
            return self._list()
        raise FoamParseError(f"{self.source}: expected list or dictionary")

    def _dict(self) -> FoamDict:
        self._skip()
        if self.text[self.i] != "{":
            raise FoamParseError(f"{self.source}: expected '{{'")
        self.i += 1
        data = FoamDict()
        while True:
            self._skip()
            if self.i >= self.n:
                raise FoamParseError(f"{self.source}: unclosed dictionary")
            if self.text[self.i] == "}":
                self.i += 1
                break
            self.parse_entry(data)
        return data

    def _list(self) -> FoamList:
        self._skip()
        if self.text[self.i] != "(":
            raise FoamParseError(f"{self.source}: expected '('")
        self.i += 1
        items = FoamList()
        while True:
            self._skip()
            if self.i >= self.n:
                raise FoamParseError(f"{self.source}: unclosed list")
            if self.text[self.i] == ")":
                self.i += 1
                break
            if self.text[self.i] == "{":
                items.append(self._dict())
            elif self.text[self.i] == "(":
                items.append(self._list())
            else:
                items.append(self._token())
        return items

    def _token(self):
        self._skip()
        if self.i >= self.n:
            raise FoamParseError(f"{self.source}: unexpected end of file")
        if self.text[self.i] == '"':
            return self._string()
        word = self._word()
        return _coerce(word)

    def _word(self) -> str:
        self._skip()
        start = self.i
        while self.i < self.n and self.text[self.i] not in " \t\r\n;{}()\"/":
            self.i += 1
        if self.i == start:
            raise FoamParseError(
                f"{self.source}: expected a word at position {self.i}"
            )
        return self.text[start:self.i]

    def _string(self) -> str:
        if self.text[self.i] != '"':
            raise FoamParseError(f"{self.source}: expected a string")
        self.i += 1
        chars = []
        while self.i < self.n:
            char = self.text[self.i]
            if char == "\\":
                self.i += 1
                if self.i >= self.n:
                    break
                chars.append(self.text[self.i])
            elif char == '"':
                self.i += 1
                return "".join(chars)
            else:
                chars.append(char)
            self.i += 1
        raise FoamParseError(f"{self.source}: unclosed string")

    def _expect_semicolon(self) -> None:
        self._skip()
        if self.i >= self.n or self.text[self.i] != ";":
            raise FoamParseError(f"{self.source}: expected ';' at position {self.i}")
        self.i += 1

    def _optional_semicolon(self) -> None:
        self._skip()
        if self.i < self.n and self.text[self.i] == ";":
            self.i += 1

    def _skip(self) -> None:
        while self.i < self.n:
            if self.text[self.i] in " \t\r\n":
                self.i += 1
                continue
            if self.text.startswith("//", self.i):
                self.i = self.text.find("\n", self.i)
                if self.i < 0:
                    self.i = self.n
                continue
            if self.text.startswith("/*", self.i):
                end = self.text.find("*/", self.i + 2)
                if end < 0:
                    raise FoamParseError(f"{self.source}: unclosed block comment")
                self.i = end + 2
                continue
            break


def _coerce(word: str):
    if word in ("yes", "no", "true", "false", "on", "off"):
        return word
    try:
        if any(ch in word for ch in ".eE"):
            return float(word)
        return int(word)
    except ValueError:
        return word
