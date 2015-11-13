# Algorithms-for-DNA-Sequencing
Coursera MOOC Algorithms for DNA Sequencing by Ben Langmead, PhD, Jacob Pritt

##Homework1
Score 6/7

##Homework2
Score **0/6** with the committed implementation.

There is a serious flaw in the way the exercices are built, as they require you to exactly copy paste
from the video lectures to get the correct answers to the questions. Indeed, they are only interested in
making sure your code does the exact same number of character comparisons or index hits as theirs, whilst 
the code they are showing is **extremely** inefficient for no valid reasons than sloppyness.

It would be much more productive to focus on finding the actual pattern matches or results, than trying to match 
the number of character comparisons... which `homework2.py` actually does: it does find the same answers as far as
the actual pattern matching is concerned :-)

While I have not tried to optimise nor to do any kinds of clever tricks, and while `homework2.py` does fairly more 
than expected (many validations, tests, etc.) in an unoptimised way, it's still much faster than the expected code. 

Question | Answer expected from Quiz | Answer from ProgrammingHomework2
-------- | --------------------------| ---------------------------------
1        |  799954                   |    85000
2        |  984143                   |   104266
3        |  127974                   |    15291
4        |   ?                       |       2
5        |  90                       |       3
6        |  79                       |       3
RunTime  |  6.847s                   |  1.917s and 1.159s w/o tests

##Homework3
Score 4/4

##Homework4
Score 4/4
