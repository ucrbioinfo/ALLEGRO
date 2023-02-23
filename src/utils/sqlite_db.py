# Class imported by ALLEGRO. No need to run it manually.
import os
import sqlite3


class GuideSQLite3DB:
    def __init__(self, db_name: str, overwrite: bool = True) -> None:
        if overwrite and os.path.exists('../src/utils/' + db_name + '.db'):
            os.remove('../src/utils/' + db_name + '.db')

        self.db_name = db_name
        self.is_connected = False
        self.connect()


    def connect(self) -> None:
        '''
        ## Args:
            db_name: The database name to create. 
                A .db file will be created in the current directory.
            overwrite: Defaults to True. If True and if `db_name` already exists, remove the old database.
        '''
        if not self.is_connected:
            self.connection = sqlite3.connect('../src/utils/' + self.db_name + '.db')
            self.connection.row_factory = self.dict_factory
            self.cursor = self.connection.cursor()
            self.is_connected = True


    def dict_factory(self, cursor, row) -> dict:
        fields = [column[0] for column in cursor.description]
        return {key: value for key, value in zip(fields, row)}


    def save_to_db(self, table: str, dictionary: dict) -> None:
        columns = ', '.join(dictionary.keys())
        placeholders = ':' + ', :'.join(dictionary.keys())

        self.cursor.execute("CREATE TABLE IF NOT EXISTS {table}(" + columns + ")")
        query = 'INSERT INTO guides ({columns}) VALUES({placeholders})'.format(
            table=table, columns=columns, placeholders=placeholders,
            )

        self.cursor.execute(query, dictionary)
        self.connection.commit()


    def save_many_to_db(self, table: str, list_of_dicts: list[dict]) -> None:
        columns = ', '.join(list_of_dicts[0].keys())
        placeholders = ':' + ', :'.join(list_of_dicts[0].keys())

        self.cursor.execute("CREATE TABLE IF NOT EXISTS {table}(" + columns + ")")
        query = 'INSERT INTO guides ({columns}) VALUES({placeholders})'.format(
            table=table, columns=columns, placeholders=placeholders,
            )

        self.cursor.executemany(query, list_of_dicts)
        self.connection.commit()

    
    def query(self, table: str, column: str, entry) -> list[dict]:
        if not self.is_connected:
            self.connect()
        
        res = self.cursor.execute("SELECT * FROM {table} WHERE {column} = '{entry}'".format(
            table=table,
            column=column,
            entry=entry
        ))

        return res.fetchall()


    def insert_coverset(self, table: str, columns: list[str], unique_columns: list[str], seq: float, score: float, covers: int) -> None:
        columns_to_str = ', '.join(columns)
        unique_columns_to_str = ', '.join(unique_columns)
        
        self.cursor.execute("CREATE TABLE IF NOT EXISTS {table}({cols}, unique ({unique_cols}))".format(
            table=table, cols=columns_to_str, unique_cols=unique_columns_to_str,
        ))

        data = (seq, score, covers)
    
        self.cursor.execute("UPDATE {table} SET {score_column} = {new_score_val} WHERE {seq_column} = {seq}".format(
            table=table,
            score_column=columns[1],
            new_score_val=score,
            seq_column=columns[0],
            seq=seq
        ))

        self.cursor.execute("INSERT OR IGNORE INTO {table} VALUES(?, ?, ?)".format(table=table), data)
        self.connection.commit()


    def close_connection(self) -> None:
        self.connection.close()
        self.is_connected = False