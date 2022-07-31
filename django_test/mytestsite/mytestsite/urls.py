"""mytestsite URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/3.1/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf.urls import url
from temp.views import get_point_ajax,get_ajax,gene,first,first1,enter,about,human,elements
from django.contrib.staticfiles.urls import staticfiles_urlpatterns


urlpatterns = [
    url('admin/', admin.site.urls),
    url('first/',first),
    url('first1/',first1),
    url('enter/',enter),
    url('about/',about),
    url('human/',human),
    url('elements/',elements),
    url('gene/',gene),
    url('get_ajax/',get_ajax),
    url('get_point_ajax/',get_point_ajax),
    
]
urlpatterns += staticfiles_urlpatterns()