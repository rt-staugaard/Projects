import sqlite3
import re
from datetime import date

conn = sqlite3.connect('japanese_study.db')
read_cursor = conn.cursor()
write_cursor = conn.cursor()
    

write_cursor.execute("DROP TABLE IF EXISTS sentence_review")

write_cursor.execute('''
            CREATE TABLE IF NOT EXISTS sentence_review(
                japanese_key TEXT PRIMARY KEY,
                times_seen INTEGER NOT NULL,
                interval INTEGER NOT NULL,
                ease_factor REAL NOT NULL,
                due_date TEXT
            )
        ''')

read_cursor.execute("SELECT japanese_text, status from sentences")

today = date.today().isoformat()

for japanese_text, status in read_cursor:

        if status == 0:
            values = (japanese_text, 0, 0, 2.5, today)

        elif status == 1:
            values = (japanese_text, 2, 1, 2.3, today)

        else: 
            values = (japanese_text, 6, 5, 2.5, today)

        write_cursor.execute("""
            INSERT OR REPLACE INTO sentence_review
            (japanese_key, times_seen, interval, ease_factor, due_date)
            VALUES (?, ?, ?, ?, ?)
        """, values)

conn.commit()
conn.close()




