
CC      := cc
CFLAGS  := -std=c11 -O2 -Wall -Wextra -I./include
BUILD   := build
SRC     := src
BIN     := $(BUILD)/app

OBJS := $(BUILD)/cyclotomic_cosets.o $(BUILD)/mymath.o $(BUILD)/main.o

$(BIN): $(OBJS)
	$(CC) $(OBJS) -o $@ -lm

$(BUILD)/%.o: $(SRC)/%.c | $(BUILD)
	$(CC) $(CFLAGS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

.PHONY: clean
clean:
	rm -rf $(BUILD)

# ==============================================================================
# AUTOMATIZACIÓN DE TIEMPOS PARA LA MEMORIA LATEX
# ==============================================================================
TEX_DIR := ..
TEX_FILE := $(TEX_DIR)/tiempos_ejecucion.tex

.PHONY: tiempos

tiempos: $(BIN_TEXT) $(BIN_ONE) $(BIN)
	@echo "Ejecutando pruebas y exportando tiempos a LaTeX..."
	@echo "% Archivo generado automaticamente por 'make tiempos'. NO EDITAR." > $(TEX_FILE)

	@echo "Midiendo tiempos de app_text..."
	@# Ejecuta app, filtra la linea del tiempo, coge la columna 5 (el numero) y lo guarda en una variable LaTeX
	@T_TEXT=$$(./$(BIN_TEXT) | grep "TIEMPO TOTAL" | awk '{print $$5}'); \
	echo "\\newcommand{\\tiempoAppText}{$$T_TEXT}" >> $(TEX_FILE)

	@echo "Midiendo tiempos de app_one..."
	@# ATENCION: Cambiar argumentos a eleccion
	@T_ONE=$$(./$(BIN_ONE) parA.txt parB.txt | grep "TIEMPO TOTAL" | awk '{print $$5}'); \
	echo "\\newcommand{\\tiempoAppOne}{$$T_ONE}" >> $(TEX_FILE)

	@echo "Midiendo tiempos de app (optimizada)..."
	@# ATENCION: Cambiar igual que antes
	@T_BIN=$$(./$(BIN) parA.txt parB.txt | grep "TIEMPO TOTAL" | awk '{print $$5}'); \
	echo "\\newcommand{\\tiempoAppOpt}{$$T_BIN}" >> $(TEX_FILE)

	@echo "¡Tiempos actualizados con éxito en $(TEX_FILE)!"
