#!/bin/bash

API_URL="https://api.statbank.dk/v1/data"

QUERY='{"table": "EJ56",
        "format": "BULK",
        "variables": [
            {"code": "OMRÅDE","values": ["000"]},
            {"code": "EJENDOMSKATE","values":["0111"]},
            {"code": "TAL","values":["100"]},
            {"code": "TID","values":["*"]}
            ]
        }'

echo "Efterspørger data fra Danmarks Statistik..."

curl -s -X POST $API_URL -H "Content-Type: application/json" -d "$QUERY" | sed 's/;/ /g' > housing_data.csv

echo "Færdig! Data er gemt i housing_data.csv"
head -n 5 housing_data.csv
