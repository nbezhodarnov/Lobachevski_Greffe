import numpy as np

def getCountsOfDigits(number):
	return abs(str(number).find('.') - len(str(number))) - 1

def polynom(a, x):
	result = 0
	for i in range(a.size):
		result += a[i] * (x ** (a.size - i - 1))
	return result

def sign(x):
	if x > 0:
		return 1
	elif x < 0:
		return -1
	else:
		return 0

def Lobachevski_Greffe_method(a, error):
	a_temp = np.zeros(a.size, dtype = np.float128)
	for i in range(a.size):
		a_temp[i] = a[i]
	a_temp_0 = np.zeros(a.size, dtype = np.float128)
	a_next = np.zeros(a.size, dtype = np.float128)
	for i in range(a.size):
		a_next[i] = a[i]
	a_status = np.zeros(a.size, dtype = int)
	need_while = True
	count = 0
	print(count, ': ', a_next)
	while(need_while):
		for i in range(a.size):
			k = 1
			add = 0
			while((i - k >= 0) and (i + k < a.size) and (k <= i)):
				add += ((-1) ** k) * a_temp[i - k] * a_temp[i + k]
				k += 1
			add *= 2
			a_next[i] *= a_next[i]
			a_next[i] += add
			if(sign(a_next[i]) == 0):
				a_temp_0[i] = a_temp[i]
			elif(sign(a_temp[i] != 0)):
				a_temp_0[i] = 0
			if(((sign(a_next[i]) != sign(a_temp[i])) and (sign(a_next[i]) != sign(a_temp_0[i]))) and (sign(a_next[i]) != 0)):
				if(a_status[i] <= 0):
					a_status[i] -= 1
					if(a_status[i] < -3):
						a_status[i] = 3
			elif(abs(a_next[i] - a_temp[i] ** 2) <= error):
				a_status[i] = 1
			elif(abs(a_next[i] - (a_temp[i] ** 2) / 2) <= error):
				a_status[i] = 2
		for i in range(a.size):
			a_temp[i] = a_next[i]
		need_while = False
		for i in range(a.size):
			if (a_status[i] <= 0):
				need_while = True
				break
		count += 1
		print(count, ': ', a_next)
	i = 1
	roots = np.zeros(a.size - 1, dtype = complex)
	complex_first = np.array([], dtype = int)
	while(i < a.size):
		if(a_status[i] == 1):
			roots[i - 1] = abs(a_next[i] / a_next[i - 1]) ** (1 / (2 ** count))
			if (abs(polynom(a, roots[i - 1])) > error):
				roots[i - 1] *= -1
		elif(a_status[i] == 2):
			roots[i - 1] = abs(a_next[i] / a_next[i - 1]) ** (1 / (2 ** (count + 1)))
			if (abs(polynom(a, roots[i - 1])) > error):
				roots[i - 1] *= -1
			roots[i] = roots[i - 1]
			i += 1
		elif(a_status[i] == 3):
			complex_first = np.append(complex_first, i - 1)
			roots[i - 1] = abs(a_next[i + 1] / a_next[i - 1]) ** (1 / (2 ** (count + 1)))
			i += 1
		i += 1
	if(complex_first.size > 0):
		j = 1
		const = -(a[1] / a[0])
		while(j < a.size):
			if(a_status[j] == 1):
				const -= roots[j - 1]
			elif(a_status[j] == 2):
				const -= 2 * roots[j - 1]
				j += 1
			elif(a_status[j] == 3):
				j += 1
			j += 1
		for k in range(complex_first.size):
			temp = const / (2 * roots[complex_first[k]])
			roots[complex_first[k]] = complex(roots[complex_first[k]].real * temp, roots[complex_first[k]].real * ((1 - temp ** 2) ** (1 / 2)))
			roots[complex_first[k] + 1] = complex(roots[complex_first[k]].real, -roots[complex_first[k]].imag)
	print()
	return roots
	
def main():
	error = 0.5e-4
	
	a = np.array([2, -1, -4, 1], dtype=np.float128)
	result = Lobachevski_Greffe_method(a, error)
	print('Roots:')
	for i in range(result.size):
		if(result[i].imag == 0):
			print(i + 1, ': ', result[i].real)
		else:
			print(i + 1, ': ', result[i].real, ' + ', result[i].imag, ' * i')
	
	
if __name__ == '__main__':
    main()
