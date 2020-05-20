FROM rocker/tidyverse

RUN sudo apt-get update && apt-get install -y \
	python3 \
	python3-pip \
	&& install2.r --error \
		--deps TRUE \
		mgcv \
		reticulate \
		gratia \
		ggExtra
