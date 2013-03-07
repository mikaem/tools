import os

ldir = os.listdir(os.path.join(os.getcwd(), 'voluviz'))
for i, l in enumerate(ldir):
	ldir[i] = eval(l)
ldir.sort()

f = open('nektovoluviz.in', 'w')
f.write(str(len(ldir)) + '\n')
for l in ldir:
    f.write(str(l) + '\n')

f.close()

