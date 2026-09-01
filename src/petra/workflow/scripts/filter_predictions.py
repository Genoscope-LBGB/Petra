"""Keep the sequences of a PETRA prediction table that satisfy an expression.

The expression is one quoted argument holding conditions on the columns of the
table, joined by and, or, not and parentheses:

    petra-filter -i predictions.final.tsv -o kept.tsv \
        "median_coverage<40 and gc_percent>59 and not phylum==Actinopteri"
"""

import argparse
import operator
import sys
from collections.abc import Callable, Iterator, Mapping, Sequence as SequenceABC
from contextlib import contextmanager
from pathlib import Path
from typing import NamedTuple

Row = Mapping[str, str]
Predicate = Callable[[Row], bool]

_COMPARISONS: dict[str, Callable[[object, object], bool]] = {
    "==": operator.eq,
    "!=": operator.ne,
    "<": operator.lt,
    "<=": operator.le,
    ">": operator.gt,
    ">=": operator.ge,
}
# Longest first at any position, so '<=' is never read as '<' then '='.
_SYMBOLS = ("<=", ">=", "==", "!=", "<", ">", "=")
_ORDERINGS = frozenset({"<", "<=", ">", ">="})
_AND = "and"
_OR = "or"
_NOT = "not"
_OPEN = "("
_CLOSE = ")"
_CONDITION = "condition"
_KEYWORDS = frozenset({_AND, _OR, _NOT})
_QUOTES = "'\""
_LENGTH_COLUMN = "sequence_length"

_EPILOG = """\
expression:
  A condition is an attribute, an operator among == != < <= > >= (= is read as
  ==), and a value: median_coverage<40. Attribute names are the column names of
  the table. Quote a value holding spaces: species=="Chelon auratus".

  A condition whose value is a number compares numerically. A sequence whose
  column holds N/A rather than a number is then dropped by ==, <, <=, > and >=,
  and kept by !=: median_coverage<40 drops a sequence of unknown coverage,
  median_coverage!=40 keeps it.

  Any other value compares as text, N/A included. Text
  cannot be ordered, so phylum>Chordata is refused.

  Conditions combine with and, or, not (also written !) and parentheses. 

examples:
  petra-filter -i predictions.final.tsv -o kept.tsv "median_coverage<40"
  petra-filter -i predictions.final.tsv "gc_percent>59 and not phylum==Actinopteri"
  petra-filter -i predictions.final.tsv -o kept.tsv \\
      "(domain==Bacteria or domain==N/A) and sequence_length>=10000"
"""


class _Token(NamedTuple):
    kind: str
    text: str


def _as_number(text: str) -> float | None:
    try:
        return float(text)
    except ValueError:
        return None


def _read_word(expression: str, position: int) -> tuple[str, bool, int]:
    """Read one word, letting quotes hold the spaces of a value."""
    letters: list[str] = []
    quoted = False
    while position < len(expression):
        character = expression[position]
        if character.isspace() or character in (_OPEN, _CLOSE):
            break
        if character in _QUOTES:
            closing = expression.find(character, position + 1)
            if closing < 0:
                raise ValueError(f"A {character} quote is never closed.")
            letters.append(expression[position + 1 : closing])
            quoted = True
            position = closing + 1
            continue
        letters.append(character)
        position += 1
    return "".join(letters), quoted, position


def _tokenize(expression: str) -> list[_Token]:
    tokens: list[_Token] = []
    position = 0
    while position < len(expression):
        character = expression[position]
        if character.isspace():
            position += 1
        elif character in (_OPEN, _CLOSE):
            tokens.append(_Token(character, character))
            position += 1
        elif character == "!" and not expression.startswith("!=", position):
            tokens.append(_Token(_NOT, character))
            position += 1
        else:
            word, quoted, position = _read_word(expression, position)
            keyword = word.lower()
            if keyword in _KEYWORDS and not quoted:
                tokens.append(_Token(keyword, word))
            else:
                tokens.append(_Token(_CONDITION, word))
    if not tokens:
        raise ValueError("The expression holds no condition.")
    return tokens


def _split_condition(text: str) -> tuple[str, str, str]:
    """Cut a condition into its attribute, operator and value."""
    for index in range(len(text)):
        for symbol in _SYMBOLS:
            if not text.startswith(symbol, index):
                continue
            attribute = text[:index].strip()
            value = text[index + len(symbol) :].strip()
            if not attribute:
                raise ValueError(f"'{text}' names no attribute before '{symbol}'.")
            if not value:
                raise ValueError(f"'{text}' compares '{attribute}' to nothing.")
            return attribute, "==" if symbol == "=" else symbol, value
    raise ValueError(
        f"'{text}' is not a condition: expected an attribute, an operator among "
        "== != < <= > >= and a value, such as median_coverage<40."
    )


def _comparison(text: str, attribute: str, symbol: str, value: str) -> Predicate:
    compare = _COMPARISONS[symbol]
    number = _as_number(value)
    if number is None:
        if symbol in _ORDERINGS:
            raise ValueError(
                f"Condition '{text}' orders '{attribute}' against the "
                f"non-numeric value '{value}': use == or != to match text."
            )
        return lambda row: compare(row[attribute], value)

    def matches(row: Row) -> bool:
        held = _as_number(row[attribute])
        if held is None:
            return symbol == "!="
        return compare(held, number)

    return matches


def _both(left: Predicate, right: Predicate) -> Predicate:
    return lambda row: left(row) and right(row)


def _either(left: Predicate, right: Predicate) -> Predicate:
    return lambda row: left(row) or right(row)


def _negate(inner: Predicate) -> Predicate:
    return lambda row: not inner(row)


class _Expression:
    """Read a filter expression as one row predicate.

    not binds tighter than and, which binds tighter than or, and parentheses
    regroup them.
    """

    def __init__(self, expression: str) -> None:
        self._tokens = _tokenize(expression)
        self._position = 0
        self.attributes: set[str] = set()

    def parse(self) -> Predicate:
        predicate = self._parse_or()
        trailing = self._peek()
        if trailing is not None:
            if trailing.kind == _CLOSE:
                raise ValueError("A ')' has no matching '('.")
            raise ValueError(
                f"Unexpected '{trailing.text}' after a complete expression: "
                "conditions are joined by 'and' or 'or'."
            )
        return predicate

    def _peek(self) -> _Token | None:
        if self._position >= len(self._tokens):
            return None
        return self._tokens[self._position]

    def _kind(self) -> str | None:
        token = self._peek()
        return None if token is None else token.kind

    def _take(self) -> _Token:
        token = self._peek()
        if token is None:
            raise ValueError(
                f"The expression ends after '{self._tokens[-1].text}': "
                "a condition is missing."
            )
        self._position += 1
        return token

    def _parse_or(self) -> Predicate:
        predicate = self._parse_and()
        while self._kind() == _OR:
            self._take()
            predicate = _either(predicate, self._parse_and())
        return predicate

    def _parse_and(self) -> Predicate:
        predicate = self._parse_not()
        while self._kind() == _AND:
            self._take()
            predicate = _both(predicate, self._parse_not())
        return predicate

    def _parse_not(self) -> Predicate:
        if self._kind() == _NOT:
            self._take()
            return _negate(self._parse_not())
        return self._parse_condition()

    def _parse_condition(self) -> Predicate:
        token = self._take()
        if token.kind == _OPEN:
            predicate = self._parse_or()
            if self._kind() != _CLOSE:
                raise ValueError("A '(' is never closed by a ')'.")
            self._take()
            return predicate
        if token.kind == _CLOSE:
            raise ValueError("A ')' has no matching '('.")
        if token.kind in (_AND, _OR):
            raise ValueError(f"'{token.text}' has no condition on its left.")
        attribute, symbol, value = _split_condition(token.text)
        self.attributes.add(attribute)
        return _comparison(token.text, attribute, symbol, value)


@contextmanager
def _destination(destination: Path | None) -> Iterator:
    if destination is None:
        yield sys.stdout
        return
    with destination.open("w", encoding="utf-8") as handle:
        yield handle


def _percentage(part: int, whole: int) -> str:
    if whole == 0:
        return "0.0%"
    return f"{part / whole:.1%}"


def filter_predictions(source: Path, destination: Path | None, expression: str) -> None:
    """Write the rows of a prediction table that satisfy an expression."""
    conditions = _Expression(expression)
    predicate = conditions.parse()

    total = retained = total_bases = retained_bases = 0
    with source.open(encoding="utf-8") as input_handle:
        header_line = input_handle.readline()
        if not header_line.strip():
            raise ValueError(
                f"{source} holds no header: expected a first line of column names."
            )
        header = header_line.rstrip("\r\n").split("\t")
        unknown = sorted(conditions.attributes.difference(header))
        if unknown:
            raise ValueError(
                f"Unknown attribute(s) {', '.join(unknown)}: the table holds "
                f"{', '.join(header)}."
            )
        measures_bases = _LENGTH_COLUMN in header

        with _destination(destination) as output_handle:
            output_handle.write("\t".join(header) + "\n")
            for line_number, line in enumerate(input_handle, start=2):
                if not line.strip():
                    continue
                values = line.rstrip("\r\n").split("\t")
                if len(values) != len(header):
                    raise ValueError(
                        f"Invalid sequence at line {line_number}: expected "
                        f"{len(header)} columns, found {len(values)}."
                    )
                row = dict(zip(header, values))
                length = (
                    int(_as_number(row[_LENGTH_COLUMN]) or 0) if measures_bases else 0
                )
                total += 1
                total_bases += length
                if predicate(row):
                    retained += 1
                    retained_bases += length
                    output_handle.write("\t".join(values) + "\n")

    report = f"Retained {retained} of {total} sequences ({_percentage(retained, total)})"
    if measures_bases:
        report += (
            f" and {retained_bases:,} of {total_bases:,} bases "
            f"({_percentage(retained_bases, total_bases)})"
        )
    print(f"{report}; expression: {expression}", file=sys.stderr)


def _parse_arguments(argv: SequenceABC[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__.splitlines()[0],
        epilog=_EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="TSV",
        help="Prediction table to filter, such as run-1/predictions.final.tsv.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="TSV",
        help="Table to write; standard output when omitted.",
    )
    parser.add_argument(
        "expression",
        nargs="+",
        metavar="EXPRESSION",
        help="Conditions to satisfy, as one quoted argument.",
    )
    return parser.parse_args(argv)


def main(argv: SequenceABC[str] | None = None) -> None:
    arguments = _parse_arguments(argv)
    filter_predictions(
        arguments.input, arguments.output, " ".join(arguments.expression)
    )


def run() -> None:
    """Installed entry point: report a rejected request without a traceback."""
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(f"Error: {error}")


if __name__ == "__main__":
    run()
