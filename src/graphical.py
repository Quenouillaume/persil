from numpy import inf
import matplotlib.pyplot as plt




def get_min_max(intervals):
	p_min = inf
	p_max = -inf
	for (x,y) in intervals:
		p_min = min(p_min,x)
		if y != inf:
			p_max = max(p_max,y)
	return p_min,p_max

def persistence_diagram(intervals):

	p_min, p_max = get_min_max(intervals)
	dp = p_max - p_min

	lower_limit = p_min - dp/5
	upper_limit = p_max + dp/5


	fig, ax = plt.subplots()
	ax.set_xlim(left = lower_limit,right = upper_limit)
	ax.set_ylim(bottom = lower_limit,top = upper_limit)
	ax.plot([lower_limit, upper_limit],[lower_limit , upper_limit],'b-')

	for (x,y) in intervals:
		if y == inf:
			ax.plot([x,x],[x,upper_limit],'g-')
		else:
			ax.plot([x],[y],'ro')
	plt.show()


def barcode(intervals):
	p_min, p_max = get_min_max(intervals)
	dp = p_max - p_min

	lower_limit = p_min - dp/5
	upper_limit = p_max + dp/5

	fig, ax = plt.subplots()
	ax.set_xlim(left = lower_limit,right = upper_limit)

	n = len(intervals)
	for i in range(n):
		(x,y) = intervals[i]
		if x == y:
			ax.plot([x],[y],'ro')
		elif y == inf:
			ax.plot([x,upper_limit],[i,i],'r-')
		else:
			ax.plot([x,y],[i,i],'r-')
	plt.show()
