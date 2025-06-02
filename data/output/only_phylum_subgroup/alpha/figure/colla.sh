#!/bin/bash

# Преобразование PDF в PNG
for file in *.pdf; do
  convert -density 300 "$file" "${file%.pdf}.png"
done

# Создание коллажа
montage *.png -tile 4x5 -geometry 400x400+10+10 collage.png

# Удаление всех PNG-файлов, кроме collage.png
for file in *.png; do
  if [[ "$file" != "collage.png" ]]; then
    rm "$file"
  fi
done

echo "Коллаж 'collage.png' создан, временные PNG-файлы удалены."
