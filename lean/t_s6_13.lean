import Poset13Step3

def s : Bool := step3All 6 26624 28672
theorem s_true : s = true := by native_decide
#print axioms s_true
