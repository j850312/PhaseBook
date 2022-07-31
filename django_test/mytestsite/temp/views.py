from django.shortcuts import render
import pandas as pd
import json
import numpy as np
import itertools
import matplotlib.pyplot as plt
import math
import re
import sklearn
from scipy.optimize import leastsq
from sklearn.metrics import r2_score
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt

def first(request):
    List = ['自強學堂', '渲染Json到模板']
    return render(request, 'first.html', {'List': List})

def first1(request):
    List = ['自強學堂', '渲染Json到模板']
    return render(request, 'tp1.html', {'List': List})
def about(request):
    List = ['自強學堂', '渲染Json到模板']
    return render(request, 'about.html', {'List': List})

        
def human(request):
    df1=pd.read_csv(r"C:/Users/CyLab/Desktop/result0.csv")
    df2 = pd.read_csv(r"C:/Users/CyLab/Desktop/result1.csv")
    df3=pd.read_csv(r"C:/Users/CyLab/Desktop/result2.csv")
    df4 = pd.read_csv(r"C:/Users/CyLab/Desktop/result3.csv")
    df5=pd.read_csv(r"C:/Users/CyLab/Desktop/result_471.csv")
    df6 = pd.read_csv(r"C:/Users/CyLab/Desktop/result_472.csv")
    df7=pd.read_csv(r"C:/Users/CyLab/Desktop/result_473.csv")
    df8 = pd.read_csv(r"C:/Users/CyLab/Desktop/result_474.csv")
    df_expression_47=pd.read_csv(r"C:/Users/CyLab\Desktop/Expression_BA47.txt",sep="\t")
    df_expression_11=pd.read_csv(r"C:/Users/CyLab\Desktop/Expression.txt",sep="\t")
    df_Annotation=pd.read_csv(r"C:/Users/CyLab\Desktop\Annotation.csv")
    df_donors=pd.read_csv(r"C:/Users/CyLab\Desktop\donors_info.txt",sep="\t")
    
# data preprocessing
# gene expresssion
    df_expression_47["Probe_id"] = df_Annotation["Probe Set ID"]
    df_expression_47 = df_expression_47.set_index("Probe_id")
    df_expression_11["Probe_id"] = df_Annotation["Probe Set ID"]
    df_expression_11 = df_expression_11.set_index("Probe_id")
    df_Annotation = df_Annotation.set_index("Probe Set ID")
    df_Annotation["Gene Symbol"] = df_Annotation["Gene Symbol"].map(lambda x : x.strip())

    list_gene_symbol = df_Annotation[(df_Annotation["Cytoband"].\
                                            map(lambda x : re.search(r"[\dXYM]+",x)).notna())].\
                                            sort_values(by = ["Cytoband"])["Gene Symbol"].\
                                            value_counts().keys().tolist()

# get target gene probe id which is the Maximum expression
    list_probe_id = []
    list_probe_id1 = []
    for item in list_gene_symbol:
        list_probe_id.append(get_max_expression_id(df_expression_annotation = df_Annotation,
                              df_expression = df_expression_11,
                              str_gene_name = item))
        list_probe_id1.append(get_max_expression_id(df_expression_annotation = df_Annotation,
                              df_expression = df_expression_47,
                              str_gene_name = item))                              

# sorted by Cytoband
    list_probe_id = df_Annotation.loc[list_probe_id].\
                sort_values(by = ["Cytoband"]).index.tolist()
    list_probe_id1 = df_Annotation.loc[list_probe_id1].\
                sort_values(by = ["Cytoband"]).index.tolist()                

# donor information
    df_donors = df_donors[pd.notna(df_donors["TOD"])]
    df_donors = df_donors.set_index("HU")
    df_donors = df_donors.sort_values("TOD")
    df_donors1= df_donors.copy()
    df_donors["BA_index"] = df_donors.index.map(lambda x : "B" + re.findall(r'(\d+)',"BA11")[0] +  ".%03d" % (x))
    df_donors1["BA_index"] = df_donors1.index.map(lambda x : "B" + re.findall(r'(\d+)',"BA47")[0] +  ".%03d" % (x))
    df_expression_11["probe_id"]=df_expression_11.index
    df_expression_47["probe_id"]=df_expression_47.index


    df_expression_11 = df_expression_11.loc[list_probe_id]
    df_expression_47 = df_expression_47.loc[list_probe_id1]
    df_donor=df_donors.to_json(orient="split")
    df_donor1=df_donors1.to_json(orient="split")
    df_expression_11=df_expression_11.to_json(orient="split")
    df_expression_47=df_expression_47.to_json(orient="split")
    gene=df1["gene"].tolist()
    df11 = df1.to_json(orient="split")
    df12 = df2.to_json(orient="split")
    df13 = df3.to_json(orient="split")
    df14 = df4.to_json(orient="split")
    df15 = df5.to_json(orient="split")
    df16 = df6.to_json(orient="split")
    df17 = df7.to_json(orient="split")
    df18 = df8.to_json(orient="split")
    
    return render(request,"human.html",{
    'data':json.dumps(df11),
    'data1':json.dumps(df12),
    'data2':json.dumps(df13),
    'data3':json.dumps(df14),
    'data4':json.dumps(df15),
    'data5':json.dumps(df16),
    'data6':json.dumps(df17),
    'data7':json.dumps(df18),
    'gene':json.dumps(gene),
    'df_donors':json.dumps(df_donor),
    'df_donors1':json.dumps(df_donor1),
    'df_expression_47':json.dumps(df_expression_47),
    'df_expression_11':json.dumps(df_expression_11)
    })


def elements(request):
    df1=pd.read_csv(r"C:/Users/CyLab/Desktop/result0.csv")
    gene=df1["gene"].tolist()
    df_donors=pd.read_csv(r"C:\Users\CyLab\Desktop\donors_info.txt",sep="\t")

    # donor information
    df_donors = df_donors[pd.notna(df_donors["TOD"])]
    df_donors = df_donors.set_index("HU")
    df_donors = df_donors.sort_values("TOD")

    df1 = df1.to_json(orient="split")
    df_donor=df_donors.to_json(orient="split")
    #return render_to_response('elements.html', locals())
    return render(request,"elements.html",{
    'df_donors':json.dumps(df_donor),
    'data':json.dumps(df1),
    'gene':json.dumps(gene),
    })
def enter(request):
    df1=pd.read_csv(r"C:/Users/CyLab/Desktop/data.csv")
    df_point=pd.read_csv(r"C:/Users/CyLab/Desktop/baboon_point.csv")
    df2 = pd.read_csv(r"C:/Users/CyLab/Desktop/data2.csv")
    df3 = df1.drop(columns="Unnamed: 0")
    df3 = df3.drop(columns="Unnamed: 1")
    df = pd.concat([df2, df3], axis=1)
    organization=[]#組織
    ogdic={}#取決於點選的組織的index
    gene_id=[]
    org=[]
    orgdic={}
    gene_name=[]
    genedic={}
    gndic={}

    for i in range(2,len(df.columns),6):
        organization.append(df.columns[i])
        ogdic[df.columns[i]]=i
    #for i in range(2,len(df_point.columns),12):
    #    if i>254:#資料少了1筆KIC之後
    #        i=i-1
    #    tp=re.split("\.|\_",df_point.columns[i])
    #    orgdic[tp[0]]=i
    for i in range(1,len(df["Unnamed: 0"])):
        gene_id.append(df["Unnamed: 0"][i])
        gene_name.append(df["Unnamed: 1"][i])
        genedic[df["Unnamed: 0"][i]]=i-1
        if df["Unnamed: 1"][i] in gndic:
            ans = []
            ans = gndic.get(df["Unnamed: 1"][i])
            ans.append(df["Unnamed: 0"][i])
            gndic[df["Unnamed: 1"][i]] = ans
        else:
            temp = []
            temp.append(df["Unnamed: 0"][i])
            gndic[df["Unnamed: 1"][i]] = temp

    gene_id.sort()
    gene_name=sorted(list(set(gene_name)))
    organization.sort()
    df1 = df.to_json(orient="split")
    df_point = df_point.to_json(orient="split")
    return render(request,"enter.html",{
    'genedic': json.dumps(genedic),
    'gndic': json.dumps(gndic),
    'gene_name':json.dumps(gene_name),
    'gene_id': json.dumps(gene_id),
    'organization':json.dumps(organization),
    'ogdic': json.dumps(ogdic),
    #'orgdic': json.dumps(orgdic),
    'data':json.dumps(df1)
    })

def gene(request):
    df1=pd.read_csv(r"C:/Users/CyLab/Desktop/data.csv")
    df2 = pd.read_csv(r"C:/Users/CyLab/Desktop/data2.csv")
    df3 = df1.drop(columns="Unnamed: 0")
    df3 = df3.drop(columns="Unnamed: 1")
    df = pd.concat([df2, df3], axis=1)
    organization=[]#組織
    ogdic={}#取決於點選的組織的index
    gene_id=[]
    gene_name=[]
    genedic={}
    gndic={}
    data=[]

    for i in range(2,len(df.columns),6):
        organization.append(df.columns[i])
        ogdic[df.columns[i]]=i
    for i in range(1,len(df["Unnamed: 0"])):
        gene_id.append(df["Unnamed: 0"][i])
        gene_name.append(df["Unnamed: 1"][i])
        genedic[df["Unnamed: 0"][i]]=i-1
        if df["Unnamed: 1"][i] in gndic:
            ans = []
            ans = gndic.get(df["Unnamed: 1"][i])
            ans.append(df["Unnamed: 0"][i])
            gndic[df["Unnamed: 1"][i]] = ans
        else:
            temp = []
            temp.append(df["Unnamed: 0"][i])
            gndic[df["Unnamed: 1"][i]] = temp

    gene_id.sort()
    gene_name=sorted(list(set(gene_name)))
    organization.sort()
    df1 = df.to_json(orient="split")
    return render(request,"gene.html",{
    'genedic': json.dumps(genedic),
    'gndic': json.dumps(gndic),
    'gene_name':json.dumps(gene_name),
    'gene_id': json.dumps(gene_id),
    'organization':json.dumps(organization),
    'ogdic': json.dumps(ogdic),
    'genedic': json.dumps(genedic),
    'data':json.dumps(df1)
    })



def sin_curve(amplitude, frequence, baseline, shift, phi):
    """"
    sine_curve : give a sine curve
    """
    return amplitude * math.sin(2 * math.pi * frequence * phi + shift) + baseline
@csrf_exempt
def get_point_ajax(request):
    gene=json.loads(request.POST.get("gene",None))
    tis=json.loads(request.POST.get("tis",None))
    df_point=pd.read_csv(r"C:/Users/CyLab/Desktop/baboon_point.csv")
    data_ans=[]
    
    for ge in gene:
        for org in tis:
            gene_index=df_point.index[df_point["EnsemblID"]==ge].tolist()
            #list_gene=df_point.loc[df_point["EnsemblID"]==ge]
            #data_ans.append(org)
            ans=""
            count=0
            for i in range(0,len(df_point.columns)):
                org_tp=re.split('\.|\_',df_point.columns[i])[0]
                
                if org_tp==org:
                    temp=df_point.ix[gene_index[0],i]
                    #temp=list_gene[list_gene.columns[i]]
                    if count!=11:
                        ans+=str(round(temp,4))+","
                        count+=1
                    else:
                        ans+=str(round(temp,4))
                    
            data_ans.append(ans)
                
    return  HttpResponse(json.dumps({"ans":data_ans}))
    #return   HttpResponse(json.dumps({"gene":gene}))


@csrf_exempt #human
def get_ajax(request):
    gene=request.POST.get("gene",None)
    age=json.loads(request.POST.get("age",None))
    sex=request.POST.get("sex",None)
    ba_select=request.POST.get("BA_select",None)
    list_probe_id=pd.read_csv(r"C:\Users\CyLab\Desktop\human_probe_id.csv",sep="\n")
    list_probe_id=list_probe_id["7932938"].tolist()
    df_Annotation=pd.read_csv(r"C:\Users\CyLab\Desktop\Annotation.csv")
    df_donors=pd.read_csv(r"C:\Users\CyLab\Desktop\donors_info.txt",sep="\t")
    # donor information
    df_donors = df_donors[pd.notna(df_donors["TOD"])]
    df_donors = df_donors.set_index("HU")
    df_donors = df_donors.sort_values("TOD")
    temp=sex.split(",")
    if ba_select=="BA_47":
        df_expression=pd.read_csv(r"C:\Users\CyLab\Desktop\Expression_BA47.txt",sep="\t")
        df_donors["BA_index"] = df_donors.index.map(lambda x : "B" + re.findall(r'(\d+)',"BA47")[0] +  ".%03d" % (x))
    else:
        df_expression=pd.read_csv(r"C:\Users\CyLab\Desktop\Expression.txt",sep="\t")
        df_donors["BA_index"] = df_donors.index.map(lambda x : "B" + re.findall(r'(\d+)',"BA11")[0] +  ".%03d" % (x))


    df_expression["Probe_id"] = df_Annotation["Probe Set ID"]
    df_expression = df_expression.set_index("Probe_id")
    df_Annotation = df_Annotation.set_index("Probe Set ID")
    df_Annotation["Gene Symbol"] = df_Annotation["Gene Symbol"].map(lambda x : x.strip())


# get target genes expression
    df_expression = df_expression.loc[list_probe_id]

    gene_tp=df_Annotation.loc[df_Annotation["Gene Symbol"]== gene]
    probe_id=gene_tp.index[0]    
    if int(age)==-1 and sex=="X":
        t1=df_donors
        point_y1=df_expression.loc[probe_id,df_donors["BA_index"]]           
    elif int(age)!=-1 and sex=="X":
        t1=df_donors.loc[df_donors["Age"].between(0,int(age))]
        t2=df_donors.loc[df_donors["Age"].between(int(age)+1,100)]
        point_y1=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(0,int(age)))]["BA_index"]]
        point_y2=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(int(age)+1,100))]["BA_index"]]
    if sex!="X" and int(age)==-1:
        if sex=="M" or temp[0]=="M":
            t3=df_donors.loc[df_donors["Sex"]=="M"]
            point_y3=df_expression.loc[probe_id,df_donors[(df_donors["Sex"]=="M")]["BA_index"]]
        if sex=="F" or temp[1]=="F":
            t4=df_donors.loc[df_donors["Sex"]=="F"]
            point_y4=df_expression.loc[probe_id,df_donors[(df_donors["Sex"]=="F")]["BA_index"]]
    elif sex=="M" and int(age)!=-1:
        r1=df_donors.loc[df_donors["Sex"]=="M"]
        t1=r1.loc[r1["Age"].between(0,int(age))]#男生0-範圍內
        t2=r1.loc[r1["Age"].between(int(age)+1,100)]#男生範圍-100內
        point_y1=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(0,int(age))) & (df_donors["Sex"]=="M")]["BA_index"]]
        point_y2=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(int(age)+1,100)) & (df_donors["Sex"]=="M")]["BA_index"]]
    elif  sex=="F" and int(age)!=-1:
        r2=df_donors.loc[df_donors["Sex"]=="F"]
        t3=r2.loc[r2["Age"].between(0,int(age))]#女生0-範圍內
        t4=r2.loc[r2["Age"].between(int(age)+1,100)]#女生範圍-100內
        point_y3=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(0,int(age))) & (df_donors["Sex"]=="F")]["BA_index"]]
        point_y4=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(int(age)+1,100)) & (df_donors["Sex"]=="F")]["BA_index"]]
    elif sex!="X" and int(age)!=-1:
        r1=df_donors.loc[df_donors["Sex"]=="M"]
        r2=df_donors.loc[df_donors["Sex"]=="F"]
        t1=r1.loc[r1["Age"].between(0,int(age))]#男生0-範圍內
        t2=r1.loc[r1["Age"].between(int(age)+1,100)]#男生範圍-100內
        t3=r2.loc[r2["Age"].between(0,int(age))]#女生0-範圍內
        t4=r2.loc[r2["Age"].between(int(age)+1,100)]#女生範圍-100內
        point_y1=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(0,int(age))) & (df_donors["Sex"]=="M")]["BA_index"]]
        point_y2=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(int(age)+1,100)) & (df_donors["Sex"]=="M")]["BA_index"]]
        point_y3=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(0,int(age))) & (df_donors["Sex"]=="F")]["BA_index"]]
        point_y4=df_expression.loc[probe_id,df_donors[(df_donors["Age"].between(int(age)+1,100)) & (df_donors["Sex"]=="F")]["BA_index"]]
    dic_last={}
    if 't1' in locals():
        dict_parameter_t1=sine_curve(t1,df_expression,df_Annotation,probe_id)
        dict_parameter_t1[probe_id]["point_x"]=t1["TOD"].tolist()
        dict_parameter_t1[probe_id]["point_y"]=point_y1.tolist()
        dic_last.setdefault('t1',dict_parameter_t1[probe_id])
        

    if 't2' in locals():
        dict_parameter_t2=sine_curve(t2,df_expression,df_Annotation,probe_id)
        dict_parameter_t2[probe_id]["point_x"]=t2["TOD"].tolist()
        dict_parameter_t2[probe_id]["point_y"]=point_y2.tolist()
        dic_last.setdefault('t2',dict_parameter_t2[probe_id])
    if 't3' in locals():
        dict_parameter_t3=sine_curve(t3,df_expression,df_Annotation,probe_id)
        dict_parameter_t3[probe_id]["point_x"]=t3["TOD"].tolist()
        dict_parameter_t3[probe_id]["point_y"]=point_y3.tolist()
        dic_last.setdefault('t3',dict_parameter_t3[probe_id])
    if 't4' in locals():
        dict_parameter_t4=sine_curve(t4,df_expression,df_Annotation,probe_id)
        dict_parameter_t4[probe_id]["point_x"]=t4["TOD"].tolist()
        dict_parameter_t4[probe_id]["point_y"]=point_y4.tolist()
        dic_last.setdefault('t4',dict_parameter_t4[probe_id])
    return  HttpResponse(json.dumps(dic_last,cls=NpEncoder)) 


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def sine_curve(df_donors, df_expression, df_expression_profile, 
                        probe_id):

    dict_parameter = {}


    Q1 = df_expression.loc[probe_id,df_donors["BA_index"]].quantile(.1)
    Q3 = df_expression.loc[probe_id,df_donors["BA_index"]].quantile(.3)

    amplitude = Q3 - Q1
    frequence = Q3 / Q1
    baseline = df_expression.loc[probe_id, df_donors["BA_index"]].mean()

    shift = 2 * math.pi * (df_donors[df_donors["BA_index"] == df_expression.loc[probe_id, df_donors["BA_index"]].sort_values().keys()[math.ceil(df_donors.shape[0]/2)]]["TOD"].tolist()[0])
    
    optimize_func = lambda x: [sin_curve(x[0],x[1],x[3],x[2],i/24)for i in df_donors["TOD"]]\
        - df_expression.loc[probe_id, df_donors["BA_index"]]
    parameters, test = leastsq(optimize_func, [amplitude, frequence, shift, baseline])

    dict_temp ={}
    dict_temp["gene"] = df_expression_profile.loc[probe_id,'Gene Symbol']
    dict_temp["Probe Set ID"]=probe_id
    dict_temp["amplitude"] = round(parameters[0],4)
    dict_temp["frequence"] = round(parameters[1],4)
    dict_temp["shift"] = round(parameters[2],4)
    dict_temp["baseline"] = round(parameters[3],4)

    dict_parameter[probe_id] = dict_temp

    del parameters, dict_temp

    # R2 

    final_r2 = sklearn.metrics.r2_score(
        df_expression.loc[probe_id,df_donors["BA_index"]],
                            [sin_curve(dict_parameter[probe_id]["amplitude"],
                                    dict_parameter[probe_id]["frequence"],
                                    dict_parameter[probe_id]["baseline"],
                                    dict_parameter[probe_id]["shift"],
                                    i/24)for i in df_donors["TOD"]])

    #             dict_parameter[probe_id]["guess_R2_score"] = guess_r2
    dict_parameter[probe_id]["finial_R2_score"] = round(final_r2,4)
    return dict_parameter

def get_max_expression_id(df_expression_annotation,
                          df_expression,
                          str_gene_name
                         ):
    """
    get probe id which is the maximum gene expression
    """
    return df_expression.loc[\
        df_expression_annotation[df_expression_annotation["Gene Symbol"] == str_gene_name].index].\
        mean(axis=1).\
        idxmax()
    #amplitude*sin(2拍(i+phase)*period)+base