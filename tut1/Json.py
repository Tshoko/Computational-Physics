# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 20:01:58 2020

@author: Boitshoko
"""

import json
import No_numpy


array=[5,10,20,100,1000]
if __name__ == '__main__':


    with open("tut1.json","r") as f:
        results = json.load(f)

    print(results)

    for q,r in results.items():
        #print(q)

        for ss in r:
            s = int(ss)
            
            #print("  ", ss)
            M=No_numpy.Make(s)[1]
            x=No_numpy.Make(s)[0]
            if(q=="xtx"):
                result =No_numpy.Dot(x,x)
                results[q][ss] = result
            elif q=="xtMx":
                result =No_numpy.XMX(M,x)
                results[q][ss] = result
    print(results)

    with open("tut1.json","w") as f:
        json.dump(results,f)
#print(Multiply(Make(5)[1]))
        import json
'''
+
+if __name__ == '__main__':
+
+
+    with open("tut1/data/sample-solutions.json","r") as f:
+        results = json.load(f)
+
+    print(results)
+
+    for q,r in results.items():
+        print(q)
+
+        for ss in r:
+            s = int(ss)
+            print("  ", s)
+
+            if results[q][ss] is None:
+                result = 42
+                results[q][ss] = result
+
+    print(results)
+
+    with open("tut1.json","w") as f:
+        json.dump(results,f)
'''