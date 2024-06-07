import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, confusion_matrix, ConfusionMatrixDisplay
from imblearn.over_sampling import SMOTE

df = pd.read_csv("dataset/GSE244266.csv",index_col=0)
print(df.head)
#print(df.describe())

print(df['stage'].value_counts())

df_label = df[['stage']]

plt.figure(figsize=(10, 6))
sns.histplot(data=df_label, x='stage')
plt.title('Distribution of Labels')
plt.xlabel('Class')
plt.ylabel('Frequency')
plt.show()

df_label.loc[df_label["stage"] == 2, "stage"] = 1 #Stage 2
df_label.loc[df_label["stage"] == 3, "stage"] = 2 #Stage 3 or 4

sm = SMOTE(sampling_strategy = 0.3, k_neighbors = 5, random_state = 42)
df, df_label = sm.fit_resample(df, df_label)

plt.figure(figsize=(10, 6))
sns.histplot(data=df_label, x='stage')
plt.title('Distribution of Labels')
plt.xlabel('Class')
plt.ylabel('Frequency')
plt.show()

X_train, X_test, y_train, y_test = train_test_split(df, df_label, test_size=0.2, random_state=42, stratify = df_label)

print("X_train:", X_train.shape, "X_test", X_test.shape)
print("y_train:", y_train.shape, "y_test", y_test.shape)

knnCls = KNeighborsClassifier(n_neighbors=5)
knnCls.fit(X_train, y_train)
y_pred = knnCls.predict(X_test)

print("Accuracy Score:", accuracy_score(y_test, y_pred))
print("Precision Score:", precision_score(y_test, y_pred, average='binary'))
print("Recall Score:", recall_score(y_test, y_pred, average='binary'))

cm = confusion_matrix(y_test, y_pred, labels=knnCls.classes_)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=knnCls.classes_)
disp.plot()
plt.show()