def bounce_seq(d):
	temp = []
	l = d.bounce_path().touch_points()
	for j in range(len(l)):
		if j == 0:
			temp += l[j]*[0]
		else:
			temp += (l[j]-l[j-1])*[j]
	return temp

def bounce_code(d):
	l = d.bounce_path().touch_points()
	n = l[-1]
	return [n-l[j] for j in range(len(l)-1)[::1]]