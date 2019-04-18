#!/bin/sh

pdflatex vsdp_workflow.tex 1> /dev/null
pdf2svg vsdp_workflow.pdf vsdp_workflow.svg
rm vsdp_workflow.pdf vsdp_workflow.aux vsdp_workflow.log
