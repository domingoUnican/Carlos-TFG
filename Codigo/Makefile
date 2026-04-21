
CC      := cc
CFLAGS  := -std=c11 -O2 -Wall -Wextra -I./include
BUILD   := build
SRC     := src
BIN     := $(BUILD)/app
BIN_TEXT := $(BUILD)/app_text
BIN_ONE  := $(BUILD)/app_one

COMMON_OBJS := $(BUILD)/pairs_reader.o $(BUILD)/mymath.o $(BUILD)/cyclotomic_cosets.o
OBJS      := $(COMMON_OBJS) $(BUILD)/main.o
TEXT_OBJS := $(COMMON_OBJS) $(BUILD)/main_txt.o
ONE_OBJS  := $(COMMON_OBJS) $(BUILD)/main-1.o

$(BIN): $(OBJS)
	$(CC) $(OBJS) -o $@ -lm

$(BIN_TEXT): $(TEXT_OBJS)
	$(CC) $(TEXT_OBJS) -o $@ -lm

$(BIN_ONE): $(ONE_OBJS)
	$(CC) $(ONE_OBJS) -o $@ -lm

$(BUILD)/%.o: $(SRC)/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

.PHONY: text one clean

text: $(BIN_TEXT)

one: $(BIN_ONE)

clean:
	rm -rf $(BUILD)

# ==============================================================================
# AUTOMATIZACIÓN DE TIEMPOS PARA LA MEMORIA LATEX (Usando /usr/bin/time)
# ==============================================================================
TEX_DIR := ..
TEX_FILE := $(TEX_DIR)/tiempos_ejecucion.tex

.PHONY: tiempos

tiempos: $(BIN_TEXT) $(BIN_ONE) $(BIN)
	@echo "Ejecutando pruebas y exportando tiempos a LaTeX..."
	@echo "% Archivo generado automaticamente por 'make tiempos'. NO EDITAR." > $(TEX_FILE)
	
	@echo "Midiendo tiempos de app_text (fuerza bruta)..."
	@T_TEXT=$$(/usr/bin/time -f "%e" ./$(BIN_TEXT) 2>&1 1>/dev/null); \
	printf "\\\\newcommand{\\\\tiempoAppText}{%s}\n" "$$T_TEXT" >> $(TEX_FILE)
	
	@echo "Midiendo tiempos de app_one (DFS básico)..."
	@T_ONE=$$(/usr/bin/time -f "%e" ./$(BIN_ONE) Tests/75/25.txt Tests/75/15.txt 3 5 2>&1 1>/dev/null); \
	printf "\\\\newcommand{\\\\tiempoAppOne}{%s}\n" "$$T_ONE" >> $(TEX_FILE)
	
	@echo "Midiendo tiempos de app (DFS Bitsets optimizado)..."
	@T_BIN=$$(/usr/bin/time -f "%e" ./$(BIN) Tests/75/25.txt Tests/75/15.txt 3 5 2>&1 1>/dev/null); \
	printf "\\\\newcommand{\\\\tiempoAppOpt}{%s}\n" "$$T_BIN" >> $(TEX_FILE)
	
	@echo "¡Tiempos actualizados con éxito en $(TEX_FILE)!"
