"""
Test the suffix-based naming feature for the complicator package.
"""

import pandas as pd
from complicator import complicate

# Read the example input file
in_name = "mgoh_low_high"
df_input = pd.read_csv(in_name + ".csv")

print("=" * 80)
print("Suffix-Based Naming Feature Test")
print("=" * 80)
print("\nInput CSV file (mgoh_low_high.csv):")
print(df_input[["Metal", "Ligand", "BETA_1", "BETA_2", "name_suffix_1", "name_suffix_2"]])
print()

# Run the complicate function
print("Running complicate function...")
df_out, equations, warnings, duplicates = complicate(
    df_in=df_input,
    water_model="SUPCRT92",
    print_warnings=True
)

print("\n" + "=" * 80)
print("Results:")
print("=" * 80)
print(f"\nGenerated {len(df_out)} complexes:")
print("\nComplex details:")
for idx, row in df_out.iterrows():
    print(f"\n  Complex {idx+1}:")
    print(f"    name:    {row['name']}")
    print(f"    abbrv:   {row['abbrv']}")
    print(f"    formula: {row['formula']}")

print("\n" + "=" * 80)
print("Verification:")
print("=" * 80)

# Verify that name has suffix but abbrv and formula don't
all_correct = True
for idx, row in df_out.iterrows():
    suffix = "_low" if "low" in row['name'] else ("_mid" if "mid" in row['name'] else "_high")
    expected_name_has_suffix = suffix in row['name']
    abbrv_no_suffix = suffix not in row['abbrv']
    formula_no_suffix = suffix not in row['formula']

    if expected_name_has_suffix and abbrv_no_suffix and formula_no_suffix:
        print(f"✓ Complex {idx+1}: name has suffix, abbrv and formula don't")
    else:
        print(f"✗ Complex {idx+1}: FAILED")
        all_correct = False

if all_correct:
    print("\n✓ All complexes have correct naming!")
else:
    print("\n✗ Some complexes have incorrect naming!")

print("\n" + "=" * 80)
print("Summary:")
print("=" * 80)
print("""
The suffix-based naming feature allows you to add suffixes to complex names
while keeping abbrv and formula as auto-generated names.

Columns in the input CSV:
  - name_suffix_1: Suffix for the first complex's name
  - name_suffix_2: Suffix for the second complex's name
  - name_suffix_3: Suffix for the third complex's name
  - name_suffix_4: Suffix for the fourth complex's name

Result columns:
  - name: Auto-generated name + suffix (e.g., "Mg(OH)+_low")
  - abbrv: Auto-generated name without suffix (e.g., "Mg(OH)+")
  - formula: Auto-generated formula without suffix (e.g., "Mg(OH)+")
""")
