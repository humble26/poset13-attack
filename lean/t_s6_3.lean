import Poset13Step3

def s : Bool := step3All 6 6144 8192
theorem s_true : s = true := by native_decide
#print axioms s_true
