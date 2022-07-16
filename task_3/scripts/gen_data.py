import sys
import random


N = 10


def create_test_data(count: int = 10) -> None:
	"""Create 'count' test files"""

	for i in range(count):
		filaname = 'input_' + str(i) + '.txt'
		f = open(filaname, 'w')
		m = 500 + i * 200

		for i in range(m):
			f.write(str(random.randint(0, 1000)) + '\n')

		f.close()


if __name__ == '__main__':
	create_test_data(N)
	sys.exit(0)
