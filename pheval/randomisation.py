import pandas as pd
import random
import logging as log

random.seed(10)


def rand(
    self: pd.DataFrame, min_num: int, max_num: int, scrambe_factor: float
) -> float:
    """
    Numeric scrambling
    Args:
        self (pd.DataFrame): dataframe records
        min_num (int): min value from this records
        max_num (int): max value from this records
        scrambe_factor (float): scramble factor scalar

    Returns:
        float: randomized number
    """
    try:
        return self + (random.uniform(min_num, max_num) * scrambe_factor)
    except Exception as err:
        log.error(self, exc_info=err)
        return self


def ssp_randomisation(df, scramble_factor=0.5) -> pd.DataFrame:
    """
      Takes as input a pandas dataframe with the semantic similarity profile and a value between 0 and 1 (scramble_factor: 0 means no scrambling and 1 means  complete randomisation). It then randomises the above scores with a degree of the scramble_factor and returns a scrambles pandas dataframe.
    Args:
          df (pd.DataFrame): original records as dataframe came from phenotype database
          scramble_factor (float) scalar scramble factor

    Returns:
        pd.Dataframe: scrambled dataframe
    """
    columns = ["SIMJ", "IC", "SCORE"]
    for c in columns:
        min_num = df[c].min()
        max_num = df[c].max()
        df[c] = df[c].apply(rand, args=(min_num, max_num, scramble_factor))
        return df
