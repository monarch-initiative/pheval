import unittest
from unittest.mock import MagicMock, patch

import polars as pl
from duckdb.duckdb import DuckDBPyConnection

from pheval.analyse.benchmark_output_type import BenchmarkOutputTypeEnum
from pheval.analyse.generate_rank_comparisons import calculate_rank_changes


class TestCalculateRankChanges(unittest.TestCase):

    def setUp(self):
        self.mock_conn = MagicMock(spec=DuckDBPyConnection)

        self.sample_true_positive_cases = pl.LazyFrame(
            {
                "result_file": ["file1", "file2", "file3", "file4"],
                "gene_identifier": ["gene1", "gene2", "gene3", "gene4"],
                "gene_symbol": ["gene1", "gene2", "gene3", "gene4"],
                "run1": [0, 2, 3, 1],
                "run2": [1, 0, 1, 2],
            }
        )

        self.run_identifiers = ["run1", "run2"]

    @patch("pheval.analyse.generate_rank_comparisons.write_table")
    @patch("pheval.analyse.generate_rank_comparisons.get_logger")
    def test_calculate_rank_changes(self, mock_get_logger, mock_write_table):
        calculate_rank_changes(
            self.mock_conn,
            self.run_identifiers,
            self.sample_true_positive_cases,
            BenchmarkOutputTypeEnum.GENE.value,
        )
        mock_get_logger.return_value.info.assert_any_call("Comparing rank changes: run1 vs. run2")
        mock_write_table.assert_called_once()
        written_df = mock_write_table.call_args[0][1]
        expected_columns = [
            "result_file",
            "gene_identifier",
            "gene_symbol",
            "run1",
            "run2",
            "rank_change",
        ]
        self.assertCountEqual(written_df.columns, expected_columns)
        self.assertTrue(
            written_df.equals(
                pl.DataFrame(
                    [
                        {
                            "result_file": "file1",
                            "gene_identifier": "gene1",
                            "gene_symbol": "gene1",
                            "run1": 0,
                            "run2": 1,
                            "rank_change": "GAINED",
                        },
                        {
                            "result_file": "file2",
                            "gene_identifier": "gene2",
                            "gene_symbol": "gene2",
                            "run1": 2,
                            "run2": 0,
                            "rank_change": "LOST",
                        },
                        {
                            "result_file": "file3",
                            "gene_identifier": "gene3",
                            "gene_symbol": "gene3",
                            "run1": 3,
                            "run2": 1,
                            "rank_change": "2",
                        },
                        {
                            "result_file": "file4",
                            "gene_identifier": "gene4",
                            "gene_symbol": "gene4",
                            "run1": 1,
                            "run2": 2,
                            "rank_change": "-1",
                        },
                    ]
                )
            )
        )
