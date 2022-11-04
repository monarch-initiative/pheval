maps = {
    "HP_HP_MAPPINGS": {
        "MAPPING_ID": "INTEGER",
        "HP_ID": "VARCHAR(10)",
        "HP_TERM": "VARCHAR(200)",
        "HP_ID_HIT": "VARCHAR(10)",
        "HP_HIT_TERM": "VARCHAR(200)",
        "SIMJ": "DOUBLE",
        "IC": "DOUBLE",
        "SCORE": "DOUBLE",
        "LCS_ID": "VARCHAR(20)",
        "LCS_TERM": "VARCHAR(150)",
    },
    "HP_MP_MAPPINGS": {
        "MAPPING_ID": "INTEGER",
        "HP_ID": "VARCHAR(10)",
        "HP_TERM": "VARCHAR(200)",
        "ZP_ID": "VARCHAR(10)",
        "ZP_TERM": "VARCHAR(200)",
        "SIMJ": "DOUBLE",
        "IC": "DOUBLE",
        "SCORE": "DOUBLE",
        "LCS_ID": "VARCHAR(40)",
        "LCS_TERM": "VARCHAR(150)",
    },
    "HP_ZP_MAPPINGS": {
        "MAPPING_ID": "INTEGER",
        "HP_ID": "VARCHAR(10)",
        "HP_TERM": "VARCHAR(200)",
        "ZP_ID": "VARCHAR(10)",
        "ZP_TERM": "VARCHAR(200)",
        "SIMJ": "DOUBLE",
        "IC": "DOUBLE",
        "SCORE": "DOUBLE",
        "LCS_ID": "VARCHAR(40)",
        "LCS_TERM": "VARCHAR(150)",
    },
}


def create_mapping(table_name):
    d = {}
    d[
        table_name
    ] = f"""CREATE TABLE EXOMISER.{table_name}_SCRAMBLE ({", ".join("{!s} {!s}".format(k, v) for (k, v) in maps[table_name].items())});"""
    return "".join(list(d.values()))


def insert_mapping(table_name):
    d = {}
    d[
        table_name
    ] = f"""INSERT INTO EXOMISER.{table_name}_SCRAMBLE ({", ".join("{!s}".format(k) for (k,v) in maps[table_name].items())}) VALUES({", ".join("?".format(k) for (k) in maps[table_name].items())});"""
    return "".join(list(d.values()))
