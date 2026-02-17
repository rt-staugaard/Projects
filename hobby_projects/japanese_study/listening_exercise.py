import sqlite3
import os
import random


def study_session():
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DB_PATH = os.path.join(BASE_DIR, "database", "japanese_study.db")
    
    conn = sqlite3.connect(DB_PATH)
    cursor = conn.cursor()

    print("--- Welcome to Oral Recognition! ---")
    print("1. New Cards (Unseen)")
    print("2. Review (Learning)")
    print("3. Known Cards")


    mode = input("Chose a mode: ")
    allowed = {"1", "2", "3"} 
    while (mode not in allowed):
        mode = input("Chose a mode: ")

    deck_size = input("Select a deck size: ")
    
    if not deck_size.isnumeric():
        deck_size = 10

    if mode == "1":
        cursor.execute(f"SELECT rowid, * FROM sentences WHERE status = 0")
    elif mode == "2":
        cursor.execute(f"SELECT rowid, * FROM sentences WHERE status = 1")
    elif mode == "3":
        cursor.execute(f"SELECT rowid, * FROM sentences WHERE status = 2")

    deck = cursor.fetchall()

    if not deck:
        print("No cards found in this category!")
        return
    
    random.shuffle(deck)

    AUDIO_PATH = os.path.join(BASE_DIR, "assets", "soundfiles")
    
    for card in deck[:int(deck_size)]:
        row_id, japanese, english, audio1, audio2, _ , _ = card
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

        if (mode == "1"):
            cursor.execute("UPDATE sentences SET status = 1 WHERE rowid = ?", (row_id,))
        if (mode == "2"): 
            was_card_difficult = input("\nDid you have any difficulty with the card? (y/n): ")
            if (was_card_difficult.lower() == "n"):
                cursor.execute("UPDATE sentences SET status = 2 WHERE rowid = ?", (row_id,))

        conn.commit()

    conn.close()
    print("Good study session! じゃあね！")

if __name__ == "__main__":
    study_session()

