from .cli import main
import os
import logging as log


def setup_logger(name, log_file, level=log.INFO):
    formatter = log.Formatter(
        fmt="%(asctime)8s - %(levelname)8s - %(module)4s - %(message)s"
    )

    handler = log.FileHandler(log_file, mode="a")
    handler.setFormatter(formatter)

    logger = log.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)
    return logger


if __name__ == "__main__":
    setup_logger("info", log_file=f"{os.path.dirname(__file__)}/../info.log")
    setup_logger(
        "debug", log_file=f"{os.path.dirname(__file__)}/../debug.log", level=log.DEBUG
    )
    main()
