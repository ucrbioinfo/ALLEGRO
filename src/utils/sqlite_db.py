import os
import sqlite3


class GuideSQLite3DB:
    def __init__(self, db_name: str, overwrite: True) -> None:
        if overwrite and os.path.exists('src/utils/' + db_name + '.db'):
            os.remove('src/utils/' + db_name + '.db')

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
            self.connection = sqlite3.connect('src/utils/' + self.db_name + '.db')
            self.connection.row_factory = self.dict_factory
            self.cursor = self.connection.cursor()
            self.is_connected = True


    def dict_factory(self, cursor, row) -> dict:
        fields = [column[0] for column in cursor.description]
        return {key: value for key, value in zip(fields, row)}


    def save_to_db(self, dictionary: dict) -> None:
        columns = ', '.join(dictionary.keys())
        placeholders = ':' + ', :'.join(dictionary.keys())

        self.cursor.execute("CREATE TABLE IF NOT EXISTS guides(" + columns + ")")
        query = 'INSERT INTO guides ({columns}) VALUES({placeholders})'.format(
            columns=columns, placeholders=placeholders,
            )

        self.cursor.execute(query, dictionary)
        self.connection.commit()


    def save_many_to_db(self, list_of_dicts: list[dict]) -> None:
        columns = ', '.join(list_of_dicts[0].keys())
        placeholders = ':' + ', :'.join(list_of_dicts[0].keys())

        self.cursor.execute("CREATE TABLE IF NOT EXISTS guides(" + columns + ")")
        query = 'INSERT INTO guides ({columns}) VALUES({placeholders})'.format(
            columns=columns, placeholders=placeholders,
            )

        self.cursor.executemany(query, list_of_dicts)
        self.connection.commit()

    
    def query(self, column: str, entry) -> list[dict]:
        if not self.is_connected:
            self.connect()
        
        res = self.cursor.execute("SELECT * FROM " + self.db_name + " WHERE " + column + " = '" + entry + "'")

        return res.fetchall()


    def close_connection(self) -> None:
        self.connection.close()
        self.is_connected = False