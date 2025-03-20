import logging

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] [%(levelname)s] [%(filename)s:%(lineno)d] - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def get_logger(name="PHEVAL"):
    return logging.getLogger(name)


def print_ascii_banner():
    """Prints ASCII banner only once when the script starts."""
    if not getattr(logging, "_ascii_printed", False):
        logging._ascii_printed = True
        pheval_banner = """
        Welcome to:
        ██████╗ ██╗  ██╗███████╗██╗   ██╗ █████╗ ██╗
        ██╔══██╗██║  ██║██╔════╝██║   ██║██╔══██╗██║
        ██████╔╝███████║█████╗  ██║   ██║███████║██║
        ██╔═══╝ ██╔══██║██╔══╝  ╚██╗ ██╔╝██╔══██║██║
        ██║     ██║  ██║███████╗ ╚████╔╝ ██║  ██║███████╗
        ╚═╝     ╚═╝  ╚═╝╚══════╝  ╚═══╝  ╚═╝  ╚═╝╚══════╝
        A framework for the empirical evaluation of phenotype-driven prioritisation tools.
        """
        print(pheval_banner)


def initialise_context(ctx):
    ctx.ensure_object(dict)
    if not getattr(ctx, "ascii_printed", False):
        ctx.ascii_printed = True
        print_ascii_banner()
