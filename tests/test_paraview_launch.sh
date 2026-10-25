#!/bin/sh
# Stub the Windows tools used by ./bin/rdesim-wsl open.
# This does not start a real paraview.exe.
set -eu

ROOT="$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)"
TIME="$ROOT/cases/annulus2D/9e-06"
if [ -e "$TIME" ]; then
    echo "refusing to replace existing $TIME" >&2
    exit 1
fi

WORK="$(mktemp -d)"
STUB_BIN="$WORK/bin"
STUB_LOG="$WORK/powershell.log"
OUT="$WORK/open.out"
EXE="$WORK/Program Files/ParaView 6.2.0/bin/paraview.exe"
mkdir -p "$STUB_BIN" "$(dirname "$EXE")"
: >"$EXE"
: >"$STUB_LOG"

cleanup() {
    rm -rf "$TIME" "$WORK"
}
trap cleanup EXIT
mkdir "$TIME"

cat >"$STUB_BIN/wslpath" <<'EOF'
#!/bin/sh
flag="$1"
path="$2"
case "$flag" in
    -w)
        case "$path" in
            *paraview.exe)
                printf '%s\n' 'C:\Program Files\ParaView 6.2.0\bin\paraview.exe'
                ;;
            *)
                rel="${path#/}"
                rel="$(printf '%s' "$rel" | sed 's|/|\\|g')"
                printf '%s\n' "\\\\wsl.localhost\\Ubuntu-26.04\\${rel}"
                ;;
        esac
        ;;
    -u)
        printf '%s\n' "$path"
        ;;
    *)
        echo "wslpath stub: unexpected arguments: $*" >&2
        exit 1
        ;;
esac
EOF

cat >"$STUB_BIN/powershell.exe" <<'EOF'
#!/bin/sh
printf '%s\n' "$*" >>"${STUB_LOG:?}"
exit 0
EOF

cat >"$STUB_BIN/tasklist.exe" <<'EOF'
#!/bin/sh
if [ "${STUB_PARAVIEW_RUNNING:-0}" = "1" ]; then
    printf '%s\n' "paraview.exe                    1234 Console                    1     10,000 K"
else
    printf '%s\n' "INFO: No tasks are running which match the specified criteria."
fi
exit 0
EOF

chmod 755 "$STUB_BIN/wslpath" "$STUB_BIN/powershell.exe" "$STUB_BIN/tasklist.exe"

fail() {
    echo "FAIL: $*" >&2
    echo "---- script output ----" >&2
    cat "$OUT" >&2 || true
    exit 1
}

require() {
    if ! grep -F -q -- "$1" "$OUT"; then
        fail "missing [$1]"
    fi
}

reject() {
    if grep -F -q -- "$1" "$OUT"; then
        fail "unexpected [$1]"
    fi
}

run_open() {
    : >"$OUT"
    : >"$STUB_LOG"
    set +e
    env -u RDESIM_PARAVIEW \
        WSL_DISTRO_NAME=Ubuntu-26.04 \
        RDESIM_PARAVIEW_WAIT=0 \
        RDESIM_PARAVIEW_EXE="$RDESIM_PARAVIEW_EXE" \
        STUB_LOG="$STUB_LOG" \
        STUB_PARAVIEW_RUNNING="${STUB_PARAVIEW_RUNNING:-0}" \
        PATH="$STUB_BIN:/usr/bin:/bin" \
        sh "$ROOT/bin/rdesim-wsl" open >"$OUT" 2>&1
    code=$?
    set -e
}

decode_first() {
    encoded="$(awk '
        {
            for (i = 1; i <= NF; i++) {
                if ($i == "-EncodedCommand" && (i + 1) <= NF) {
                    print $(i + 1)
                    exit
                }
            }
        }
    ' "$STUB_LOG")"
    [ -n "$encoded" ] || return 1
    printf '%s' "$encoded" | base64 -d | iconv -f UTF-16LE -t UTF-8
}

RDESIM_PARAVIEW_EXE="$EXE"
STUB_PARAVIEW_RUNNING=1
run_open
[ "$code" -eq 0 ] || fail "tasklist success exited $code"
require "Windows reports paraview.exe is running."
require "Exit code:"
require "ParaView executable: $EXE"
reject "ParaView is running (pid"
decoded="$(decode_first)" || fail "could not decode Start-Process"
printf '%s\n' "$decoded" | grep -F -q "Program Files" || fail "decoded command dropped Program Files"
printf '%s\n' "$decoded" | grep -F -q "ParaView 6.2.0" || fail "decoded command dropped ParaView 6.2.0"
printf '%s\n' "$decoded" | grep -F -q "annulus2D.foam" || fail "decoded command dropped annulus2D.foam"
printf '%s\n' "$decoded" | grep -F -q "Ubuntu-26.04" || fail "decoded command dropped Ubuntu-26.04"
printf '%s\n' "$decoded" | grep -F -q "Start-Process -FilePath" || fail "decoded command is not Start-Process"
echo "ok: tasklist confirms paraview.exe"

STUB_PARAVIEW_RUNNING=0
run_open
[ "$code" -ne 0 ] || fail "missing paraview.exe was treated as success"
require "Windows Task Manager has no paraview.exe"
require "Command:"
require "Exit code:"
require 'If paraview.exe starts and no window appears, reinstall ParaView on Windows.'
require '\\wsl.localhost\Ubuntu-26.04\'
require '\\wsl$\Ubuntu-26.04\'
require "annulus2D.foam"
reject "Windows reports paraview.exe is running."
echo "ok: tasklist miss prints the command and the paste path"

RDESIM_PARAVIEW_EXE="/no/such/paraview.exe"
STUB_PARAVIEW_RUNNING=1
run_open
[ "$code" -ne 0 ] || fail "a missing RDESIM_PARAVIEW_EXE exited 0"
require "RDESIM_PARAVIEW_EXE is set but that file was not found: /no/such/paraview.exe"
reject "Windows reports paraview.exe is running."
if [ -s "$STUB_LOG" ]; then
    fail "a missing executable still started PowerShell"
fi
echo "ok: missing RDESIM_PARAVIEW_EXE fails without a launch"
