import pandas as pd
import math
from complicator import complicate

# Test case 1: BETA_1 and BETA_3 present, but BETA_2 missing
# Expected: Only BETA_1 should be estimated, BETA_3 should be skipped
test_data_1 = {
    "Metal": ["Mg+2"],
    "Ligand": ["Cl-"],
    "BETA_1": [1.5],
    "BETA_2": [float('NaN')],
    "BETA_3": [3.0],
    "BETA_4": [float('NaN')],
}

# Test case 2: BETA_1, BETA_2, and BETA_4 present, but BETA_3 missing
# Expected: BETA_1 and BETA_2 should be estimated, BETA_4 should be skipped
test_data_2 = {
    "Metal": ["Mg+2"],
    "Ligand": ["Cl-"],
    "BETA_1": [1.5],
    "BETA_2": [2.5],
    "BETA_3": [float('NaN')],
    "BETA_4": [4.0],
}

# Test case 3: BETA_2 and BETA_3 present, but BETA_1 missing
# Expected: Nothing should be estimated (BETA_2 and BETA_3 both require BETA_1)
test_data_3 = {
    "Metal": ["Mg+2"],
    "Ligand": ["Cl-"],
    "BETA_1": [float('NaN')],
    "BETA_2": [2.5],
    "BETA_3": [3.0],
    "BETA_4": [float('NaN')],
}

print("=" * 80)
print("Test Case 1: BETA_1 and BETA_3 present, BETA_2 missing")
print("Expected: Only first complex estimated, third complex skipped")
print("=" * 80)
df_test_1 = pd.DataFrame(test_data_1)
df_out_1, equations_1, warnings_1, duplicates_1 = complicate(
    df_in=df_test_1,
    water_model="SUPCRT92",
    print_warnings=True
)
print(f"\nNumber of complexes estimated: {len(df_out_1)}")
print(f"Warnings: {warnings_1}")
print()

print("=" * 80)
print("Test Case 2: BETA_1, BETA_2, BETA_4 present, BETA_3 missing")
print("Expected: First and second complexes estimated, fourth complex skipped")
print("=" * 80)
df_test_2 = pd.DataFrame(test_data_2)
df_out_2, equations_2, warnings_2, duplicates_2 = complicate(
    df_in=df_test_2,
    water_model="SUPCRT92",
    print_warnings=True
)
print(f"\nNumber of complexes estimated: {len(df_out_2)}")
print(f"Warnings: {warnings_2}")
print()

print("=" * 80)
print("Test Case 3: BETA_2 and BETA_3 present, BETA_1 missing")
print("Expected: No complexes estimated (all require BETA_1)")
print("=" * 80)
df_test_3 = pd.DataFrame(test_data_3)
df_out_3, equations_3, warnings_3, duplicates_3 = complicate(
    df_in=df_test_3,
    water_model="SUPCRT92",
    print_warnings=True
)
print(f"\nNumber of complexes estimated: {len(df_out_3)}")
print(f"Warnings: {warnings_3}")
print()

print("=" * 80)
print("All tests completed!")
print("=" * 80)
