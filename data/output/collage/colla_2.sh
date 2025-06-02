#!/bin/bash

# Проверка, передан ли аргумент (путь к папке)
if [ -z "$1" ]; then
  echo "Использование: $0 <путь_к_папке>"
  exit 1
fi

TARGET_DIR="$1"

# Переход в указанную папку
cd "$TARGET_DIR" || {
  echo "Не удалось перейти в папку: $TARGET_DIR"
  exit 1
}

# Создание временной папки для коллажа
TEMP_DIR="collage_temp"
mkdir -p "$TEMP_DIR"

# Преобразование PDF в PNG
for file in *.pdf; do
  convert -density 300 "$file" "${file%.pdf}.png"
done

# Создание коллажа во временной папке
montage *.png -tile 4x8 -geometry 400x400+10+10 "$TEMP_DIR/collage.png"

# Проверка успешности создания
if [ -f "$TEMP_DIR/collage.png" ]; then
  # Удаление всех PNG-файлов, кроме коллажа
  rm -f *.png
  # Перемещение коллажа обратно
  mv "$TEMP_DIR/collage.png" .
  # Удаление временной папки
  rmdir "$TEMP_DIR"
  echo "Коллаж 'collage.png' создан и сохранён в: $TARGET_DIR"
else
  echo "Ошибка: collage.png не был создан."
  exit 2
fi

