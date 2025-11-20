import pandas as pd
import os
from complicator import complicate

# Create test data with custom names for Mg-OH complexes
test_data = {
    "Metal": ["Mg+2", "Mg+2", "Mg+2"],
    "Ligand": ["OH-", "OH-", "OH-"],
    "BETA_1": [2.56, 2.56, 2.56],
    "BETA_2": [5.02, 5.02, 5.02],
    "BETA_3": [float('NaN'), float('NaN'), float('NaN')],
    "BETA_4": [float('NaN'), float('NaN'), float('NaN')],
    "name_1": ["MgOH+_low", "MgOH+_mid", "MgOH+_high"],
    "name_2": ["Mg(OH)2_low", "Mg(OH)2_mid", "Mg(OH)2_high"],
}

df_test = pd.DataFrame(test_data)

print("=" * 80)
print("Testing custom naming feature")
print("=" * 80)
print("\nInput dataframe:")
print(df_test[["Metal", "Ligand", "BETA_1", "BETA_2", "name_1", "name_2"]])
print()

# Run complicate for each row
for idx in range(len(df_test)):
    print(f"\n--- Processing row {idx+1} ---")
    df_single = df_test.iloc[[idx]]

    df_out, equations, warnings, duplicates = complicate(
        df_in=df_single,
        water_model="SUPCRT92",
        print_warnings=False
    )

    if not df_out.empty:
        print(f"Generated complex names:")
        for i, name in enumerate(df_out["name"].values):
            print(f"  Complex {i+1}: {name}")
    else:
        print("No complexes generated")

print("\n" + "=" * 80)
print("Test completed!")
print("=" * 80)

# Also test that auto-naming still works when custom names are not provided
print("\n" + "=" * 80)
print("Testing auto-naming (no custom names provided)")
print("=" * 80)

test_data_auto = {
    "Metal": ["Ca+2"],
    "Ligand": ["Cl-"],
    "BETA_1": [0.5],
    "BETA_2": [0.8],
}

df_test_auto = pd.DataFrame(test_data_auto)
print("\nInput dataframe:")
print(df_test_auto)
print()

df_out_auto, _, _, _ = complicate(
    df_in=df_test_auto,
    water_model="SUPCRT92",
    print_warnings=False
)

if not df_out_auto.empty:
    print(f"Generated complex names (auto):")
    for i, name in enumerate(df_out_auto["name"].values):
        print(f"  Complex {i+1}: {name}")

print("\n" + "=" * 80)
print("All tests completed!")
print("=" * 80)
