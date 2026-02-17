import json
import sqlite3
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DB_PATH = os.path.join(BASE_DIR,"database","japanese_study.db")
JSON_PATH = os.path.join(BASE_DIR,"data_raw","kanji_with_readings.json")
SIMILARITY_PATH = os.path.join(BASE_DIR,"data_raw","kanji_similarity_list.txt")


def import_kanji():

    with open(JSON_PATH,'r',encoding='utf-8') as f:
        data1 = json.load(f)

        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()

        cursor.execute("DROP TABLE IF EXISTS kanji_meanings")

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS kanji_meanings(
                kanji TEXT PRIMARY KEY,
                meanings TEXT,
                onyomi TEXT,
                kunyomi TEXT,
                stroke_count INTEGER, 
                jlpt INTEGER
            )
        ''')

        kanji_data = data1.get("kanjis", {})

        kanji_list = []
        for kanji_char, info in kanji_data.items():
            kanji_list.append((
                kanji_char, 
                ", ".join(info.get('meanings', [])), 
                ", ".join(info.get('on_readings', [])), 
                ", ".join(info.get('kun_readings', [])), 
                info.get('stroke_count'), 
                info.get('jlpt')
            ))
            
        cursor.executemany("INSERT OR IGNORE INTO kanji_meanings VALUES (?, ?, ?, ?, ?, ?)", kanji_list)

        conn.commit()
        print(f"Successfully imported {len(kanji_list)} kanji!")

    with open(SIMILARITY_PATH,'r',encoding='utf-8') as f:
        data2 = json.load(f)

        conn = sqlite3.connect(DB_PATH)
        cursor = conn.cursor()

        cursor.execute("DROP TABLE IF EXISTS kanji_similarities")

        cursor.execute('''
            CREATE TABLE IF NOT EXISTS kanji_similarities(
                kanji TEXT PRIMARY KEY,
                similarity_group INTEGER,
                group_size INTEGER
            )
        ''')

        kanji_list = []
        group_nr = 0
        for group in data2:
            group_nr += 1
            group_size = len(group)

            for kanji in group:
                kanji_list.append((kanji,group_nr,group_size    
            ))
            
        cursor.executemany("INSERT OR IGNORE INTO kanji_similarities VALUES (?, ?, ?)", kanji_list)

        conn.commit()
        print(f"Successfully imported {len(kanji_list)} kanji!")
    
    cursor.execute("DROP TABLE IF EXISTS full_kanji_table")

    cursor.execute('''
        CREATE TABLE full_kanji_table AS
        SELECT 
            m.kanji, m.meanings, m.onyomi, m.kunyomi, m.stroke_count, m.jlpt,
            s.similarity_group, s.group_size
        FROM kanji_meanings m
        RIGHT JOIN kanji_similarities s ON m.kanji = s.kanji
        ''')

    conn.commit()
    print(f"Complete the combined table!")
    conn.close()

if __name__ == "__main__":
    import_kanji()
            
            



