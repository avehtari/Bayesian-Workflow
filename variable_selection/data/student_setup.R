
library(rprojroot)
root<-has_file(".ARMM-Examples-root")$make_fix_file()

# read mathematics and Portuguese language results
d1=read.table(root("Student/data","student-mat.csv"),sep=";",header=TRUE)
d2=read.table(root("Student/data","student-por.csv"),sep=";",header=TRUE)

# merge mathematics and Portuguese language results
d3=merge(d1,d2,by=c("school","sex","age","address","famsize","Pstatus","Medu","Fedu","Mjob","Fjob","reason","nursery","internet"),all=TRUE)
print(nrow(d3)) # 382 students

# pick unique columns
data=d3[,c("G1.x","G2.x","G3.x","G1.y","G2.y","G3.y","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime.y","studytime.y","failures.y","schoolsup.y","famsup.y","paid.y","activities.y", "nursery", "higher.y", "internet", "romantic.y","famrel.y","freetime.y","goout.y","Dalc.y","Walc.y","health.y","absences.y")]

# rename columns
colnames(data)<-c("G1mat","G2mat","G3mat","G1por","G2por","G3por","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime","studytime","failures","schoolsup","famsup","paid","activities", "nursery", "higher", "internet", "romantic","famrel","freetime","goout","Dalc","Walc","health","absences")

# for simplicity transform binary factors to 0/1 coding
predictors_binary <- c("school","sex","address","famsize","Pstatus","schoolsup","famsup","paid","activities","nursery","higher","internet","romantic")
data[,predictors_binary] <- apply(t(1:13), 2, function (x) { as.numeric(data[,predictors_binary[x]])-1 })

# save merged data
write.csv(data, root("Student/data","student-merged-all.csv"),row.names=FALSE)

# merge mathematics and Portuguese language results
d3=merge(d1,d2,by=c("school","sex","age","address","famsize","Pstatus","Medu","Fedu","Mjob","Fjob","reason","nursery","internet"),all=TRUE)
print(nrow(d3)) # 682 students

# pick unique columns
data=d3[,c("G1.x","G2.x","G3.x","G1.y","G2.y","G3.y","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime.y","studytime.y","failures.y","schoolsup.y","famsup.y","paid.y","activities.y", "nursery", "higher.y", "internet", "romantic.y","famrel.y","freetime.y","goout.y","Dalc.y","Walc.y","health.y","absences.y")]

# rename columns
colnames(data)<-c("G1mat","G2mat","G3mat","G1por","G2por","G3por","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime","studytime","failures","schoolsup","famsup","paid","activities", "nursery", "higher", "internet", "romantic","famrel","freetime","goout","Dalc","Walc","health","absences")

# for simplicity transform binary factors to 0/1 coding
predictors_binary <- c("school","sex","address","famsize","Pstatus","schoolsup","famsup","paid","activities","nursery","higher","internet","romantic")
data[,predictors_binary] <- apply(t(1:13), 2, function (x) { as.numeric(data[,predictors_binary[x]])-1 })

# save merged data
write.csv(data, root("Student/data","student-merged-all.csv"),row.names=FALSE)

datapor=d2[,c("G1","G2","G3","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime","studytime","failures","schoolsup","famsup","paid","activities", "nursery", "higher", "internet", "romantic","famrel","freetime","goout","Dalc","Walc","health","absences")]

# rename columns
colnames(datapor)<-c("G1por","G2por","G3por","school","sex","age","address","famsize","Pstatus","Medu","Fedu","traveltime","studytime","failures","schoolsup","famsup","paid","activities", "nursery", "higher", "internet", "romantic","famrel","freetime","goout","Dalc","Walc","health","absences")

# for simplicity transform binary factors to 0/1 coding
predictors_binary <- c("school","sex","address","famsize","Pstatus","schoolsup","famsup","paid","activities","nursery","higher","internet","romantic")
datapor[,predictors_binary] <- apply(t(1:13), 2, function (x) { as.numeric(datapor[,predictors_binary[x]])-1 })

# save merged data
write.csv(datapor, root("Student/data","student-por-numeric.csv"),row.names=FALSE)
