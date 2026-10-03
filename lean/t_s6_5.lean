import Poset13Step3

def s : Bool := step3All 6 10240 12288
theorem s_true : s = true := by native_decide
#print axioms s_true
