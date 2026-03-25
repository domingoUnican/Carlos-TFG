
CC      := cc
CFLAGS  := -std=c11 -g -O0 -Wall -Wextra -I./include
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
