#!/usr/bin/env python3
"""Strip C/C++/CUDA comments from a source file; write the result to stdout.

Used by EncodeKernelFiles.cmake so kernel doc/comment text does not get baked
into the CudaKernelSources / CommonKernelSources string literals. Those literals
ship in the .so and are handed to NVRTC on every cold JIT (and are part of the
SHA1 the JIT cache is keyed on), so comments are pure runtime payload there.

Line numbers are preserved: every newline inside a removed block comment is kept,
and line comments keep their terminating newline. This keeps NVRTC diagnostics,
ncu's source view, and any source cross-reference aligned with the .cu/.cc file.

A small state machine (not a regex) so that '//' or '/*' appearing inside string
or character literals -- e.g. a URL or a printf format -- is left untouched.
Known limitation: C++ raw string literals (R"(...)") are not modeled; the kernel
sources do not use them.
"""

import sys


def strip(src: str) -> str:
    out = []
    i, n = 0, len(src)
    state = "code"  # code | line | block | string | char
    while i < n:
        c = src[i]
        d = src[i + 1] if i + 1 < n else ""
        if state == "code":
            if c == "/" and d == "/":
                # A comment becomes one space (C standard) so it still separates
                # tokens: a//x\nb -> "a \nb", never "ab".
                state = "line"
                out.append(" ")
                i += 2
            elif c == "/" and d == "*":
                state = "block"
                out.append(" ")
                i += 2
            elif c == '"':
                state = "string"
                out.append(c)
                i += 1
            elif c == "'":
                state = "char"
                out.append(c)
                i += 1
            else:
                out.append(c)
                i += 1
        elif state == "line":
            if c == "\n":
                state = "code"
                out.append(c)  # keep the newline
            i += 1
        elif state == "block":
            if c == "*" and d == "/":
                state = "code"
                i += 2
            else:
                if c == "\n":
                    out.append(c)  # preserve line count
                i += 1
        elif state == "string":
            out.append(c)
            if c == "\\" and d:
                out.append(d)
                i += 2
            else:
                if c == '"':
                    state = "code"
                i += 1
        elif state == "char":
            out.append(c)
            if c == "\\" and d:
                out.append(d)
                i += 2
            else:
                if c == "'":
                    state = "code"
                i += 1
    return "".join(out)


def main() -> int:
    if len(sys.argv) != 2:
        sys.stderr.write("usage: strip_comments.py <file>\n")
        return 2
    with open(sys.argv[1], "r", encoding="utf-8", errors="surrogateescape") as f:
        src = f.read()
    sys.stdout.write(strip(src))
    return 0


if __name__ == "__main__":
    sys.exit(main())
