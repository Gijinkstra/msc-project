MAIN = main
AUX  = build/bin
OUT  = build

.PHONY: all clean

all:
	mkdir -p $(AUX)
	latexmk -pdf -output-directory=$(AUX) $(MAIN).tex
	cp $(AUX)/$(MAIN).pdf $(OUT)/$(MAIN).pdf

clean:
	rm -rf $(AUX)
	rm -f $(OUT)/$(MAIN).pdf