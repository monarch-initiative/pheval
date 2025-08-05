import unittest

import polars as pl

from pheval.post_processing.validate_result_format import ResultSchema, validate_dataframe


class TestResultSchema(unittest.TestCase):
    def setUp(self):
        """Set up valid test DataFrames for each schema."""
        self.valid_gene_df = pl.DataFrame(
            {
                "gene_symbol": ["BRCA1", "TP53"],
                "gene_identifier": ["ENSG00000012048", "ENSG00000141510"],
                "score": [0.9, 0.8],
                "grouping_id": ["group1", "group2"],
            }
        )

        self.valid_variant_df = pl.DataFrame(
            {
                "chrom": ["1", "X"],
                "start": [100, 200],
                "end": [150, 250],
                "ref": ["A", "C"],
                "alt": ["G", "T"],
                "score": [0.95, 0.85],
                "grouping_id": ["group1", "group2"],
            }
        )

        self.valid_disease_df = pl.DataFrame(
            {
                "disease_name": ["Cancer", "Alzheimer's"],
                "disease_identifier": ["DOID:162", "DOID:10652"],
                "score": [0.97, 0.92],
                "grouping_id": ["group1", "group2"],
            }
        )

    def test_valid_gene_schema(self):
        """Test validation for a correct gene result schema."""
        self.assertTrue(ResultSchema.GENE_RESULT_SCHEMA.validate(self.valid_gene_df))

    def test_missing_column(self):
        """Test that missing a required column raises ValueError."""
        invalid_df = self.valid_gene_df.drop("gene_symbol")
        with self.assertRaises(ValueError, msg="Missing required column: gene_symbol"):
            ResultSchema.GENE_RESULT_SCHEMA.validate(invalid_df)

    def test_incorrect_type(self):
        """Test that a column with the wrong data type raises TypeError."""
        invalid_df = self.valid_gene_df.with_columns(pl.col("score").cast(pl.Utf8))
        with self.assertRaises(TypeError, msg="Column 'score' has type Utf8, expected Float64"):
            ResultSchema.GENE_RESULT_SCHEMA.validate(invalid_df)

    def test_null_grouping_id(self):
        """Test that null values in `grouping_id` raise ValueError."""
        invalid_df = self.valid_gene_df.with_columns(pl.lit(None).alias("grouping_id"))
        with self.assertRaises(ValueError, msg="'grouping_id' column should not contain null values if provided."):
            ResultSchema.GENE_RESULT_SCHEMA.validate(invalid_df)

    def test_missing_grouping_id_column(self):
        """Test that missing `grouping_id` does not raise an error."""
        valid_df_without_grouping = self.valid_gene_df.drop("grouping_id")
        self.assertTrue(ResultSchema.GENE_RESULT_SCHEMA.validate(valid_df_without_grouping))

    def test_valid_variant_schema(self):
        """Test validation for a correct variant result schema."""
        self.assertTrue(ResultSchema.VARIANT_RESULT_SCHEMA.validate(self.valid_variant_df))

    def test_valid_disease_schema(self):
        """Test validation for a correct disease result schema."""
        self.assertTrue(ResultSchema.DISEASE_RESULT_SCHEMA.validate(self.valid_disease_df))


class TestValidateDataframeDecorator(unittest.TestCase):
    def setUp(self):
        """Set up valid test DataFrames for decorator testing."""
        self.valid_gene_df = pl.DataFrame(
            {
                "gene_symbol": ["BRCA1"],
                "gene_identifier": ["ENSG00000012048"],
                "score": [0.9],
                "grouping_id": ["group1"],
            }
        )

        self.mock_function = validate_dataframe(ResultSchema.GENE_RESULT_SCHEMA)(
            self.mock_function.__get__(self, TestValidateDataframeDecorator)
        )

    @staticmethod
    def mock_function(self, df: pl.DataFrame):
        """A mock function that simply returns True if executed."""
        return True

    def test_decorator_valid_df(self):
        """Test that a function with a valid DataFrame executes normally."""
        self.assertTrue(self.mock_function(self.valid_gene_df))

    def test_decorator_invalid_df(self):
        """Test that a function with an invalid DataFrame raises an error."""
        invalid_df = self.valid_gene_df.drop("gene_symbol")
        with self.assertRaises(ValueError):
            self.mock_function(invalid_df)
