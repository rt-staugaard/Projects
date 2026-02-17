#! /bin/bash

SYMBOL=["APPL , NVDA"] 
API_KEY="VFV89RNABI3POOUT"

echo "Fetching data for $SYMBOL[]..."

curl -s "https://www.alphavantage.co/query?function=TIME_SERIES_DAILY&symbol=$SYMBOL&apikey=$API_KEY&datatype=csv" \
    | tail -n +2 \
    | tr -d '\r' > $OUTPUT_FILE