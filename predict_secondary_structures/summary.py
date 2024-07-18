from sklearn.metrics import confusion_matrix

with open('compare_results.txt','r') as f:
    lines=f.readlines()
    lines=[line.strip().split('\t') for line in lines]
labels=[line[-2] for line in lines]
predicted_labels=[line[-1] for line in lines]
tn, fp, fn, tp = confusion_matrix(labels, predicted_labels).ravel()
accuracy = (tp+tn)/(tp+tn+fp+fn)
sensitivity = tp/(tp+fn) 
specificity = tn/(tn+fp)
ppv = tp/(tp+fp)
npv = tn/(tn+fn)
f1score = 2*tp/(2*tp+fp+fn)
mcc=(tp*tn-fp*fn)/(((tp+fn)*(tp+fp)*(tn+fp)*(tn+fn))**0.5)
print("accuracy: {:.4f},sensitivity: {:.4f},specificity: {:.4f},ppv: {:.4f},npv: {:.4f},F1 score: {:.4f},mcc: {:.4f}".format(accuracy,sensitivity,specificity,ppv,npv,f1score,mcc))