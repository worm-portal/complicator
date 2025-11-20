"""
Demonstration of the custom naming feature for the complicator package.

This script shows how to use custom names for complexes by providing
name_1, name_2, name_3, and/or name_4 columns in your input CSV file.
"""

import os
import pandas as pd
from complicator import complicate

# Read the example input file
in_name = "mgoh_low_high"
df_input = pd.read_csv(in_name + ".csv")

print("=" * 80)
print("Custom Naming Feature Demonstration")
print("=" * 80)
print("\nInput CSV file (mgoh_low_high.csv):")
print(df_input)
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
print(f"\nGenerated {len(df_out)} complexes with custom names:")
for idx, row in df_out.iterrows():
    print(f"  - {row['name']}")

print("\n" + "=" * 80)
print("Summary:")
print("=" * 80)
print("""
The custom naming feature allows you to specify custom names for each complex
by adding columns to your input CSV file:
  - name_1: Custom name for the first complex
  - name_2: Custom name for the second complex
  - name_3: Custom name for the third complex
  - name_4: Custom name for the fourth complex

If you don't provide a custom name, the default auto-generated name will be used.

In this example, we created three sets of Mg-OH complexes with different
association constants and gave them descriptive names:
  - MgOH+_low, Mg(OH)2_low
  - MgOH+_mid, Mg(OH)2_mid
  - MgOH+_high, Mg(OH)2_high
""")
