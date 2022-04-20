import random

file1 = open("test.txt", "w")
for i in range(100):
    file1.write(str(random.randint(0, pow(10, 12))) + "\n")