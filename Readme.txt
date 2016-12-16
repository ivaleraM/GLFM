{\rtf1\ansi\ansicpg1252\cocoartf1348\cocoasubrtf170
{\fonttbl\f0\fnil\fcharset0 HelveticaNeue;\f1\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red38\green38\blue38;\red66\green66\blue66;}
\margl1440\margr1440\vieww23620\viewh8060\viewkind0
\deftab720
\pard\pardeftab720\sl640\sa320

\f0\b\fs54 \cf2 \expnd0\expndtw0\kerning0
General Table Completion using a Bayesian Nonparametric Model
\f1\fs28 \
\pard\pardeftab720\sl500\sa320
\cf2 \expnd0\expndtw0\kerning0
Compile\

\b0 \cf2 \expnd0\expndtw0\kerning0
Before using the code, please add path 'C\\ code' and compile the sampler function in Matlab with the command:\
mex  -lgsl -lgmp IBPsampler.cpp\
\pard\pardeftab720\sl500\sa320

\b \cf2 \expnd0\expndtw0\kerning0
Table completion function\

\b0 \expnd0\expndtw0\kerning0
The function matrix_completion.m completes the table with the MAP solution for the missing data\cf0 \expnd0\expndtw0\kerning0
. Please refer to the function for details on the inputs and outputs of the function.\cf2 \expnd0\expndtw0\kerning0
\
\pard\pardeftab720\sl500\sa320

\b \cf2 \expnd0\expndtw0\kerning0
Reproduce experiments in the paper \

\b0 \cf2 \expnd0\expndtw0\kerning0
The script simMiss.m reproduce the experiment on the \cf0 \expnd0\expndtw0\kerning0
QSAR biodegradation database.\cf2 \expnd0\expndtw0\kerning0
\
\pard\pardeftab720\sl500\sa320

\b \cf2 \expnd0\expndtw0\kerning0
Description\

\b0 \expnd0\expndtw0\kerning0
Please refer to Description.m\cf0 \kerning1\expnd0\expndtw0 \
\pard\pardeftab720\sl500\sa320

\b \cf2 \expnd0\expndtw0\kerning0
Citation\
\pard\pardeftab720\sl512\sa320

\b0 \cf2 \expnd0\expndtw0\kerning0
Please, cite it as detailed below.\cf0 \kerning1\expnd0\expndtw0 \
\pard\pardeftab720\sl420
\cf3 \expnd0\expndtw0\kerning0
I. Valera and Z. Ghahramani, 
\b \expnd0\expndtw0\kerning0
"General Table Completion using a Bayesian Nonparametric Model"
\b0 \expnd0\expndtw0\kerning0
, Neural Information Processing Systems Conference 2014 (NIPS 2014). Montreal (Canada), 2014.\
\
\pard\pardeftab720\sl500\sa320

\b \cf2 \expnd0\expndtw0\kerning0
Contact
\b0 \cf3 \expnd0\expndtw0\kerning0
\
\pard\pardeftab720\sl420
\cf3 \expnd0\expndtw0\kerning0
For further information or contact Isabel Valera at ivalera@mpi-sws.org.\
}