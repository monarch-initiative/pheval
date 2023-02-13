from click import Option, UsageError


class InputError(Exception):
    """Exception raised for missing required inputs."""

    def __init__(self, file, message="Missing required input"):
        self.file: str = file
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} -> {self.file} "


class MutuallyExclusiveOptionError(Option):
    """Exception raised for when"""

    def __init__(self, *args, **kwargs):
        self.mutually_exclusive = set(kwargs.pop("mutually_exclusive", []))
        help_ = kwargs.get("help", "")
        if self.mutually_exclusive:
            ex_str = ", ".join(self.mutually_exclusive)
            kwargs["help"] = help_ + (
                " NOTE: This argument is mutually exclusive with " " arguments: [" + ex_str + "]."
            )
        super(MutuallyExclusiveOptionError, self).__init__(*args, **kwargs)

    def handle_parse_result(self, ctx, opts, args):
        if self.mutually_exclusive.intersection(opts) and self.name in opts:
            raise UsageError(
                "Illegal usage: `{}` is mutually exclusive with "
                "arguments `{}`.".format(self.name, ", ".join(self.mutually_exclusive))
            )

        return super(MutuallyExclusiveOptionError, self).handle_parse_result(ctx, opts, args)


class IncorrectFileFormatError(Exception):
    def __init__(self, file, expectation, message="Incorrect File Type"):
        self.file: str = file
        self.expectation: str = expectation
        self.message: str = message
        super().__init__(self.message)

    def __str__(self):
        return f"{self.message} -> {self.file} (expected {self.expectation})"
