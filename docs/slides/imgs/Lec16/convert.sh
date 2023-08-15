for pdf in *.pdf ; do
    sips --resampleWidth 1000 -s format png --out "${pdf%%.*}.png" "$pdf"
done