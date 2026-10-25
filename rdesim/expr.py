"""Arithmetic expressions used by optimization objectives.

Identifiers, numbers, parentheses, and + - * / only. There is no call syntax.
"""

from __future__ import annotations


class ExpressionError(Exception):
    pass


def evaluate(expression: str, names: dict[str, float]) -> float:
    parser = _Expr(_tokenize(expression), expression)
    value = parser.parse_sum()
    if parser.i != len(parser.tokens):
        raise ExpressionError(f"Unexpected token in '{expression}'")
    return float(value)


def _tokenize(expression: str) -> list[str]:
    tokens = []
    i = 0
    text = expression.strip()
    while i < len(text):
        char = text[i]
        if char.isspace():
            i += 1
            continue
        if char in "+-*/()":
            tokens.append(char)
            i += 1
            continue
        if char.isdigit() or char == ".":
            start = i
            i += 1
            while i < len(text) and (text[i].isdigit() or text[i] in ".eE+-"):
                if text[i] in "+-" and text[i - 1] not in "eE":
                    break
                i += 1
            tokens.append(text[start:i])
            continue
        if char.isalpha() or char == "_":
            start = i
            i += 1
            while i < len(text) and (text[i].isalnum() or text[i] == "_"):
                i += 1
            tokens.append(text[start:i])
            continue
        raise ExpressionError(f"Invalid character {char!r} in '{expression}'")
    return tokens


class _Expr:
    def __init__(self, tokens: list[str], source: str):
        self.tokens = tokens
        self.source = source
        self.i = 0

    def parse_sum(self) -> float:
        value = self.parse_product()
        while self._peek() in ("+", "-"):
            op = self._pop()
            rhs = self.parse_product()
            value = value + rhs if op == "+" else value - rhs
        return value

    def parse_product(self) -> float:
        value = self.parse_unary()
        while self._peek() in ("*", "/"):
            op = self._pop()
            rhs = self.parse_unary()
            if op == "/":
                if rhs == 0:
                    raise ExpressionError(f"Division by zero in '{self.source}'")
                value = value / rhs
            else:
                value = value * rhs
        return value

    def parse_unary(self) -> float:
        if self._peek() in ("+", "-"):
            op = self._pop()
            value = self.parse_unary()
            return value if op == "+" else -value
        return self.parse_atom()

    def parse_atom(self) -> float:
        token = self._pop()
        if token is None:
            raise ExpressionError(f"Unexpected end of '{self.source}'")
        if token == "(":
            value = self.parse_sum()
            if self._pop() != ")":
                raise ExpressionError(f"Missing ')' in '{self.source}'")
            return value
        try:
            return float(token)
        except ValueError:
            raise ExpressionError(
                f"Unknown name '{token}' in '{self.source}'. "
                "Names must be objective entries written by the solver."
            ) from None

    def _peek(self):
        if self.i >= len(self.tokens):
            return None
        return self.tokens[self.i]

    def _pop(self):
        token = self._peek()
        if token is not None:
            self.i += 1
        return token


def evaluate_with_names(expression: str, names: dict[str, float]) -> float:
    """Evaluate after checking that every identifier is in ``names``."""
    tokens = _tokenize(expression)
    needed = [token for token in tokens if token[:1].isalpha() or token[:1] == "_"]
    missing = [name for name in needed if name not in names]
    if missing:
        raise ExpressionError(
            "Objective expression uses names that are not objective outputs: "
            + ", ".join(missing)
        )
    parser = _Expr(tokens, expression)
    # Rebind atom lookup by wrapping names. The simple evaluate() treats
    # identifiers as errors, so substitute values first.
    replaced = []
    for token in tokens:
        if token[:1].isalpha() or token[:1] == "_":
            replaced.append(format(float(names[token]), ".16g"))
        else:
            replaced.append(token)
    return _Expr(replaced, expression).parse_sum()
