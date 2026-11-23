"""The window. Everything here may import Qt; nothing outside may."""


def prepare(app) -> None:
    """Give a QApplication the look every MIL-C window is drawn and tested
    with.

    The Fusion style, on every system. Under the native macOS style the
    check boxes of a table were drawn for its first row only and form
    fields kept their narrowest width — and none of it showed in the
    tests, which run without a screen and therefore in Fusion already.
    One style everywhere means what is tested is what is seen."""
    from . import style
    from .theme import apply_defaults
    app.setStyle("Fusion")
    style.apply(app)
    apply_defaults()
