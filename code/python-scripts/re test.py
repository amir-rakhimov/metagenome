import re
test_list=['k__a',"p__aa","c__asb"]
r=re.compile("k__")
newlist=list(filter(r.match, test_list))
print(''.join(newlist))