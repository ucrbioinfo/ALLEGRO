import sqlite3

connection = sqlite3.connect("coversets.db")

def dict_factory(cursor, row):
    fields = [column[0] for column in cursor.description]
    return {key: value for key, value in zip(fields, row)}

connection.row_factory = dict_factory
cursor = connection.cursor()


# coversets = {
#     'ABCD': (23.2, {1, 2, 3}),
#     'ABED': (2.2, {4, 3})
# }

tups = [('ABCD', 3.3, 3), ('ABCD', 3.3, 3), ('ABCD', 23.3, 24)]

table_name = 'coversets'
columns = 'guide_sequence, average_score, covers'
unique_columns = 'guide_sequence, covers'

cursor.execute("CREATE TABLE IF NOT EXISTS {table}({cols}, unique ({unique_cols}))".format(
    table=table_name, cols=columns, unique_cols=unique_columns,
))

for tup in tups:
    cursor.execute("INSERT OR IGNORE INTO {table} VALUES(?, ?, ?)".format(table=table_name), tup)


cursor.execute("UPDATE {table} SET average_score = 100 WHERE guide_sequence = 'ABCD'".format(
    table=table_name
))

res = cursor.execute("SELECT * FROM {table} WHERE guide_sequence = 'ABCD'".format(
    table=table_name
))

print(res.fetchall())

connection.close()