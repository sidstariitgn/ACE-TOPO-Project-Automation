import GuiVersion
import GuiDBM
import GuiFILE
import GuiML
import GuiMO
import GuiVC
import GuiRC
import GuiSC
import GuiRun

#The following line is for backwards-compatibility, DO NOT DELETE IT.
GuiVersion.RecordVersion("2026.0.0.14417")

# =====================================================================
# SETTINGS
# =====================================================================
SOURCE_MODEL = "Hf-1(1000)"           # model with the ALD chemistry (name exactly as in the TOPO tree)
TARGET_MODEL = "Hf-HAR"         # single-material HAR model (make this the active model)

# VC numbers are model-specific: leave empty and set them by hand,
# or fill in the Oxide VC number(s) of Hf-HAR here.
OXIDE_VCS = []

SAVE_PATH = r"J:\Saketh_TOPO\3D_Nand\HAR\Hf-HAR.DTF"
SUBMIT = False                  # True -> submit after saving (set BCs first)

# Model options (same as your previous run)
END_TIME = 0.125
N_STEPS = 10000
MAX_DT = 0.3
FEATURE_TEMPERATURE = 523

# Adaption (your previous values: 0.8, 0.8, 1, 6)
REFINE_THRESHOLD = 0.9
COARSEN_RATIO = 0.85
MIN_LEVEL = 3
MAX_LEVEL = 7
# =====================================================================


def copy_item(item_type, name):
    GuiDBM.Copy(item_type, name, GuiDBM.MODEL, name, GuiDBM.MODEL,
                SOURCE_MODEL, TARGET_MODEL, True)


GuiFILE.SetMode("TOPO")

# ---- 1. Chemistry from the source model ----------------------------
copy_item(GuiDBM.SPECIES, "TDMAHF")
copy_item(GuiDBM.SPECIES, "TDMAZ")
copy_item(GuiDBM.BULKSPECIES, "HZO(B)")
for name in ["HF-INTER_1(S)", "HF-INTER_2(S)", "HZO-INTER(S)", "OH-SITE(S)",
             "TDMAHF-PHYSISORBED(S)", "TDMAZ-PHYSISORBED(S)",
             "ZR-INTER_1(S)", "ZR-INTER_2(S)"]:
    copy_item(GuiDBM.SURFACESPECIES, name)
copy_item(GuiDBM.SURFACESITE, "Surface_Speceis")
copy_item(GuiDBM.MIXTURE, "HZO")
copy_item(GuiDBM.MIXTURE, "SUB")
copy_item(GuiDBM.SURFACEREACTION, "ALD")

# ---- 2. Oxide material (same sequence as your journal:
#         create, copy to target, remove from source) ------------------
GuiDBM.AddItem(GuiDBM.MIXTURE, "New_Mixture", GuiDBM.MODEL)
GuiDBM.AddToComposition(GuiDBM.MIXTURE, "New_Mixture", "SIO2(B)", 1.000000, GuiDBM.MODEL)
GuiDBM.Rename(GuiDBM.MIXTURE, "New_Mixture", "Oxide", GuiDBM.MODEL)
copy_item(GuiDBM.MIXTURE, "Oxide")
GuiDBM.Delete(GuiDBM.MIXTURE, "Oxide", GuiDBM.MODEL, SOURCE_MODEL)

# ---- 3. Active mixtures / species, VC -------------------------------
GuiML.SetMixtures(["HZO", "Oxide"], 1)
if OXIDE_VCS:
    GuiVC.Set(OXIDE_VCS, "Chem/MixtureName", "Oxide")
GuiML.SetSpecies(["O", "O2", "TDMAHF", "TDMAZ"], 1)

# ---- 4. Model options -----------------------------------------------
GuiMO.Set("Shared/StartTime", 0)
GuiMO.Set("Shared/EndTime", END_TIME)
GuiMO.Set("Shared/NSteps", N_STEPS)
GuiMO.Set("Shared/MaxDT", MAX_DT)
GuiMO.Set("Chem/ApplyCFLLimit", 1)
GuiMO.Set("Chem/FeatureTemperature", FEATURE_TEMPERATURE)

# ---- 5. Reaction conditions -----------------------------------------
GuiRC.SetArray("Chem/ReactionMechanismName", "Oxide", "ALD")
GuiRC.SetArray("Chem/ReactionMechanismName", "HZO", "ALD")

# ---- 6. Solution control (adaption) ---------------------------------
GuiSC.Set("Relax/Threshold", REFINE_THRESHOLD)
GuiSC.Set("Relax/CoarsenRatio", COARSEN_RATIO)
GuiSC.Set("Adapt/MinLevel", MIN_LEVEL)
GuiSC.Set("Adapt/MaxLevel", MAX_LEVEL)

# ---- 7. Save / run ---------------------------------------------------
GuiFILE.SaveAs(SAVE_PATH)
if SUBMIT:
    GuiRun.Submit()