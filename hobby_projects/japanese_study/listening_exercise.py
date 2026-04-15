import sqlite3
import os
import random
from datetime import datetime, timedelta

def update_review(current_ease, current_interval, quality):
    new_ease = current_ease + (0.1 - (5 - quality) * (0.08 + (5 - quality) * 0.02))
    if new_ease < 1.3: new_ease = 1.3 

    if quality < 3: 
        new_interval = 1
    else:
        new_interval = round(current_interval * new_ease)
        
    return new_ease, new_interval

def study_session():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DB_PATH = os.path.join(BASE_DIR, "database", "japanese_study.db")
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    print("--- Welcome to Oral Recognition! ---")
    print("1. New Cards (Unseen)")
    print("2. Daily Review")


    mode = input("Chose a mode: ")
    allowed = {"1", "2"} 
    while (mode not in allowed):
        mode = input("Chose a mode: ")

    deck_size = input("Select a deck size: ")
    
    if not deck_size.isnumeric():
        deck_size = 10

    if mode == "1":
        cursor.execute(f'''SELECT s.rowid, s.japanese_text, s.english_text, s.audio_1, 
                       s.audio_2, r.ease_factor, r.interval
                       FROM sentence_review r
                       JOIN sentences s ON r.rowid = s.rowid
                       WHERE status = 0''')
    elif mode == "2":
        cursor.execute(f'''SELECT s.rowid, s.japanese_text, s.english_text, s.audio_1, 
                       s.audio_2, r.ease_factor, r.interval
                       FROM sentence_review r
                       JOIN sentences s ON r.rowid = s.rowid
                       WHERE r.due_date <= datetime('now')
                       AND status = 1
                       ORDER BY r.due_date ASC''')

    deck = cursor.fetchall()

    if not deck:
        print("No cards found in this category!")
        return

    AUDIO_PATH = os.path.join(BASE_DIR, "assets", "soundfiles")
    
    for card in deck[:int(deck_size)]:
        row_id, japanese, english, audio1, audio2, ease_factor, interval = card
        selected_audio = random.choice([a for a in [audio1, audio2] if a])

        show_answer = False
        cmd = ''
        while True:
            if (not show_answer or cmd == 'r'):
                print("\nPlaying audio...")
                os.system(f"mpv --no-video --vo=null '{AUDIO_PATH}/{selected_audio}'")

            cmd = input("\n[r] Repeat Audio | [Enter] Show Answer / Next Card | [q] Quit: ").lower().strip()
        
            if (cmd == 'r'):
                continue

            elif ((cmd == 'e' or cmd == 'edit') and show_answer == True):
                edited_translation = input(f"Displaying Current English: {english}\nNew translation: ")
                cursor.execute("UPDATE sentences SET english_text = ? WHERE rowid = ?", (edited_translation, row_id,))
                conn.commit()
                print("Update complete")
                continue

            elif ((cmd == 'j_edit') and show_answer == True):
                edited_translation = input(f"Displaying Current Japanese: {japanese}\nNew translation: ")
                cursor.execute("UPDATE sentences SET japanese_text = ? WHERE rowid = ?", (edited_translation, row_id,))
                conn.commit()
                print("Update complete")
                continue

            elif cmd == 'q':
                conn.close()
                return
            
            elif cmd == '':
                if not show_answer:
                    print(f" Japanese: {japanese}")
                    print(f" English: {english}")
                    show_answer = True
                else:
                    break
                
        difficulty = ''
        allowed = {"1", "2", "3", "4", "5"} 
        while (difficulty not in allowed):
            difficulty = input("Select difficulty: [1, 2, 3, 4, 5] (1 = Not understood, 3 = inteligble, 5 = Completely)\n")
        
        new_ease, new_interval = update_review(ease_factor, interval, int(difficulty))
        new_date_obj = datetime.now() + timedelta(days = new_interval)
        formatted_date = new_date_obj.strftime('%Y-%m-%d %H:%M:%S')

        if mode == "1":
            cursor.execute("UPDATE sentences SET status = 1 WHERE rowid = ?", (row_id,))
        
        cursor.execute('''
            UPDATE sentence_review 
            SET ease_factor = ?, interval = ?, due_date = ?
            WHERE rowid = ?
        ''', (new_ease, new_interval, formatted_date, row_id))

        conn.commit()
    conn.close()
    print("Good study session! じゃあね！")

if __name__ == "__main__":
    study_session()

