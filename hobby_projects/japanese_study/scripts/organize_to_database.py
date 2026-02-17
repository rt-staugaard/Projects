import sqlite3
import re

conn = sqlite3.connect('japanese_study.db')
cursor = conn.cursor()

cursor.execute("DROP TABLE IF EXISTS sentences")

cursor.execute('''
            CREATE TABLE IF NOT EXISTS sentences(
                japanese_text TEXT,
                english_text TEXT,
               audio_1 TEXT,
               audio_2 TEXT
            )
        ''')

audio_regex = re.compile(r"\[sound:(.*?)\]")

with open("Ankidrone Sentence Pack V4.txt","r",encoding="utf-8") as f:
    i = 0
    for line in f:
        if i < 2:
            i += 1
            continue

        parts = line.strip().split('\t')

        if len(parts) < 1:
            continue

        japanese_text = parts[0]


        english_text = ""
        for j in range(len(parts) - 1):
            current_segment = parts[j + 1]
            if current_segment == "":
                continue
            elif re.match(r"\[",current_segment):
                break
            elif "<br" in current_segment:
                clean_part = current_segment.split("<br")[0]
                english_text += clean_part
                break
            else:
                english_text += parts[j + 1]

        all_audio = audio_regex.findall(line)

        audio_1 = next((s for s in all_audio if "Male"  in s), None)
        audio_2 = next((s for s in all_audio if "Female" in s), None)
        
        if (audio_1 and audio_2 and len(all_audio) == 1):
            audio_1 = all_audio[0]

        if japanese_text and (audio_1 or audio_2):
            cursor.execute(
                "INSERT INTO sentences (japanese_text, english_text, audio_1, audio_2) VALUES (?,?,?,?)",
                (japanese_text, english_text, audio_1, audio_2) 
            )

conn.commit()
conn.close()

print("Database created successfully!")

        