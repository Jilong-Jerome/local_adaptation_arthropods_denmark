string = f"""
awk -F '\\t' '{{{{count0=0; count1=0; for(i=3; i<=NF; i++) {{{{count0+=gsub(/0/, "", $i); count1+=gsub(/1/, "", $i)}}}}; print $0, count0, count1}}}}' | \\
"""

print(string)
number = 0.01
print(f"{number}")

test_str = "#dsagdgg"
if ("##" in test_str):
    print("true")
else:
    print("false")


test_str  = "0/1/1/1/1/1/0/0"
print(test_str.count("0"))


a= [1,2,34,4,5]
print(sum(a))

test_str = f"""echo "position\\tx\\tn\\n" > spid_chro_pop.header"""
print(test_str)
