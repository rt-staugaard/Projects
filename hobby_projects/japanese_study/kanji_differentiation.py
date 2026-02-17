import sqlite3
import os
import random

class Kanji:
    def __init__(self, kanji, meaning, on_reading, kun_reading):
        self.char = kanji
        self.meaning = meaning
        self.onyomi = on_reading
        self.kunyomi = kun_reading

def study_session():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DB_PATH = os.path.join(BASE_DIR, "database", "japanese_study.db")
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    test_size = input("Select a number of test: ")
    
    if not test_size.isnumeric():
        test_size = 10

    cursor.execute("SELECT DISTINCT similarity_group FROM full_kanji_table")
    selected_entries = [row[0] for row in cursor.fetchall()]
    random.shuffle(selected_entries)

    for i in range(int(test_size)):
        cursor.execute("SELECT kanji, meanings, onyomi, kunyomi  FROM full_kanji_table WHERE similarity_group = ?", (selected_entries[i],))
        rows = cursor.fetchall()
        
        kanji_list = []
        for row in rows:
            new_kanji = Kanji(*row)
            kanji_list.append(new_kanji)

        random_entry = random.randint(0,len(kanji_list) - 1)
        selected_kanji = kanji_list[random_entry]

        while (len(kanji_list) < 4):
            readings = selected_kanji.onyomi 
            if isinstance(readings, str):
                readings = [readings]
            
            placeholders = ', '.join(['?'] * len(readings))
            query = f"SELECT kanji, meanings, onyomi, kunyomi FROM full_kanji_table WHERE onyomi IN ({placeholders})"
            cursor.execute(query, readings)
            aditional_kanji = cursor.fetchall()
            for row in aditional_kanji:
                new_kanji = Kanji(*row)
                if new_kanji.char not in [k.char for k in kanji_list]:
                    kanji_list.append(new_kanji)

            if len(kanji_list) < 4:
                cursor.execute("SELECT kanji, meanings, onyomi, kunyomi FROM full_kanji_table ORDER BY RANDOM() LIMIT 4")
                fillers = cursor.fetchall()
                for f in fillers:
                    new_kanji = Kanji(*f)
                    if new_kanji.char not in [k.char for k in kanji_list]:
                        kanji_list.append(new_kanji)
        
        kanji_list = kanji_list[0:8]
        random.shuffle(kanji_list)
        correct_index = kanji_list.index(selected_kanji) + 1

        print(f"\nWhich kanji has the following meanings: {selected_kanji.meaning}")
        n = 0
        for kanji in kanji_list:
            n += 1 
            print(n, kanji.char)

        while True:
                answer = input("\nSelect a number: ")
                if (answer.isnumeric()):
                    break    

        if (int(answer) == correct_index):
            print("Correct")
        else:
            print(f"No, the correct kanji was: {selected_kanji.char}. The selected kanji ({kanji_list[int(answer) - 1].char}) has the meanings: {kanji_list[int(answer) - 1].meaning}")

if __name__ == "__main__":
    study_session()

