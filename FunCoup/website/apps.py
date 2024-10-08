from django.apps import AppConfig

# You shouldnt have to touch anything in here
class WebsiteConfig(AppConfig):
    default_auto_field = 'django.db.models.BigAutoField'
    name = 'website'
