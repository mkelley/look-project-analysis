data := $(subst data/stacks,data/plots,$(patsubst %.fits,%.png,$(wildcard data/stacks/*/*med.fits)))

.PHONY: dirs pngs

default: $(data)

data/plots/%.png : data/stacks/%.fits
	@mkdir -p $(shell dirname $@)
	fits2bitmap --cmap=viridis --stretch=log -o $@ $<
