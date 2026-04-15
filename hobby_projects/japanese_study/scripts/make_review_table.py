import sqlite3
import re
from datetime import date

conn = sqlite3.connect('../database/japanese_study.db')
read_cursor = conn.cursor()
write_cursor = conn.cursor()
    

write_cursor.execute("DROP TABLE IF EXISTS sentence_review")

write_cursor.execute('''
    CREATE TABLE IF NOT EXISTS sentence_review(
        rowid INT PRIMARY KEY,
        interval INTEGER NOT NULL,
        ease_factor REAL NOT NULL,
        due_date TEXT NOT NULL,
        FOREIGN KEY (rowid) REFERENCES sentences(rowid)
    )
''')

write_cursor.execute("CREATE INDEX IF NOT EXISTS idx_sentence_due_date ON sentence_review(due_date)")

read_cursor.execute("SELECT rowid, status from sentences")

today = date.today().isoformat()

for rowid, status in read_cursor:

        if status == 0:
            values = (rowid, 0, 2.5, today)

        elif status == 1:
            values = (rowid, 1, 2.3, today)

        else: 
            values = (rowid, 5, 2.5, today)

        write_cursor.execute("""
            INSERT OR REPLACE INTO sentence_review
            (rowid, interval, ease_factor, due_date)
            VALUES (?, ?, ?, ?)
        """, values)

conn.commit()
conn.close()




